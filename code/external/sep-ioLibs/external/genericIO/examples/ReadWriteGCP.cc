#include <float3DReg.h>
#include <ioModes.h>
#include <stdlib.h>
#include <chrono>
using namespace SEP;
using namespace std::chrono;
int main(int argc, char **argv) {
  ioModes modes(argc, argv);
  std::shared_ptr<genericIO> io = modes.getIO("GCPBUFFERS");

  if (argc < 3) {
    std::cerr << argv[0] << "read/write  directory"
              << " [n1 n2 n3] " << std::endl;
    exit(-1);
  }
  std::string mode = argv[1];
  if (mode != std::string("read") && mode != std::string("write") &&
      mode != std::string("both")) {
    std::cerr << argv[0] << "read/write/both directory"
              << " [n1 n2 n3]" << std::endl;
    exit(-1);
  }
  std::string dir = std::string(argv[2]);
  high_resolution_clock::time_point t2, t3, t1;

  if (mode == std::string("write") || mode == std::string("both")) {
    int n1 = 1000, n2 = 1000, n3 = 200;
    if (argc >= 4) n1 = atoi(argv[3]);
    if (argc >= 5) n2 = atoi(argv[4]);
    if (argc >= 6) n3 = atoi(argv[5]);
    std::vector<int> nw(3), fw(3, 0), jw(3, 1);
    nw[0] = n1;
    nw[1] = n2;
    nw[2] = n3;

    std::shared_ptr<genericRegFile> file = io->getRegFile(dir, usageOut);

    std::shared_ptr<hypercube> hyper(new hypercube(n1, n2, n3 * 10));
    std::shared_ptr<float3DReg> buf(new float3DReg(hyper));
    for (int i = 0; i < (long long)n1 * (long long)n2 * (long long)n3; i++) {
      buf->getVals()[i] = i;
    }

    file->setHyper(hyper);
    t1 = high_resolution_clock::now();

    for (int iw = 0; iw < 10; iw++) {
      fw[2] = n3 * iw;
      file->writeFloatWindow(nw, fw, jw, buf);
    }
    file->writeDescription();

    file->close();

    t2 = high_resolution_clock::now();
    auto d2 = duration_cast<microseconds>(t2 - t1).count();
    std::cout << "To cloud " << (double)buf->getHyper()->getN123() * 4 / d2
              << " MB/s " << std::endl;
  }
  if (mode == std::string("read") || mode == std::string("both")) {
    std::shared_ptr<genericRegFile> file = io->getRegFile(dir, usageIn);
    std::shared_ptr<hypercube> hyper = file->getHyper();
    std::shared_ptr<float3DReg> buf(new float3DReg(hyper));
    t1 = high_resolution_clock::now();

    std::vector<int> nw(3, 1000), fw(3, 0), jw(3, 1);
    nw = hyper->getNs();
    file->readFloatWindow(nw, fw, jw, buf);
    file->close();
    t2 = high_resolution_clock::now();
    auto d2 = duration_cast<microseconds>(t2 - t1).count();
    std::cout << "From cloud " << (double)buf->getHyper()->getN123() * 4 / d2
              << " MB/s " << d2 << std::endl;
  }
}