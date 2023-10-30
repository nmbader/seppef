#include "functions.h"

// Functions definition

void MB::print(std::shared_ptr<float1DReg> vec, std::string message, std::ostream &o) {
    std::string espace(message.size(), ' ');
    o << message ;
    unsigned int nz = vec->getHyper()->getAxis(1).n;
    std::shared_ptr<float1D> vals = vec->_mat;

	o << "[ ";
	for (int iz =0; iz < nz; iz++){
		o << (*vals)[iz] << " ";
	}
	o << "]" << std::endl;
}

void MB::print(std::shared_ptr<float2DReg> vec, std::string message, std::ostream &o) {
    std::string espace(message.size(), ' ');
    o << message ;
    unsigned int nz = vec->getHyper()->getAxis(1).n;
    unsigned int nx = vec->getHyper()->getAxis(2).n;
    std::shared_ptr<float2D> vals = vec->_mat;
    for (int ix =0; ix < nx; ix++){
        o << "[ ";
        for (int iz =0; iz < nz; iz++){
            o << (*vals)[ix][iz] << " ";
        }
        o << "]" << std::endl;
    }
}

void MB::print(std::shared_ptr<complex2DReg> cvec, std::string message, std::ostream &o) {
    std::string espace(message.size(), ' ');
    o << message ;
    unsigned int nz = cvec->getHyper()->getAxis(1).n;
    unsigned int nx = cvec->getHyper()->getAxis(2).n;
	std::shared_ptr<complex2D> cvals = cvec->_mat;
    for (int ix =0; ix < nx; ix++){
        o << "[ ";
        for (int iz =0; iz < nz; iz++){
            o <<"("<< (*cvals)[ix][iz].real() << "+i" << (*cvals)[ix][iz].imag() <<")  ";
        }
        o << "]" << std::endl;
    }
}

std::shared_ptr<float1DReg> MB::sincWavelet(float wc, int N, float alpha){
	std::shared_ptr<float1DReg> vec (new float1DReg(2*N+1));
	std::shared_ptr<float1D> vals = vec->_mat;

	float num;
	for (int iz=0; iz<2*N+1; iz++){
		num = wc*sinc(wc*M_PI*(iz-N))*(alpha +(1-alpha)*cos(M_PI*(iz-N)/N));
		(*vals)[iz] = num;
	}
	return vec;
}

std::shared_ptr<float1DReg> MB::gaussianWavelet(float sig, int N){
	std::shared_ptr<float1DReg> vec (new float1DReg(2*N+1));
	std::shared_ptr<float1D> vals = vec->_mat;

	float num = 0;
	for (int iz=N; iz<2*N+1; iz++){ 
		num = exp(-0.5*(iz-N)*(iz-N)/(sig*sig));
		(*vals)[iz] = num;
		(*vals)[2*N-iz] = num;
	}
	return vec;
}

std::shared_ptr<float1DReg> MB::rickerWavelet(float wc, int N){
	std::shared_ptr<float1DReg> vec (new float1DReg(2*N+1));
	std::shared_ptr<float1D> vals = vec->_mat;

	float num = 0;
	float sig = sqrt(2)/(M_PI*wc);
	for (int iz=N; iz<2*N+1; iz++){ 
		num = (1-(iz-N)*(iz-N)/(sig*sig))*exp(-0.5*(iz-N)*(iz-N)/(sig*sig));
		(*vals)[iz] = num;
		(*vals)[2*N-iz] = num;
	}
	return vec;
}

void MB::readParam(int argc, char **argv, std::string param, int &value, int defaultVal){
	
	bool answer = false;
	int i = 1;
	param = param + "=";
	std::string str;

	while ( (answer == false) && (i < argc) ) {
		str = argv[i];
		if (str.substr(0, param.length()) == param){
			answer = true;
			str.erase(0, param.length());
			value = std::stoi(str);
		}
		i++;
	}

	if (answer == false)
		value = defaultVal;
}

void MB::readParam(int argc, char **argv, std::string param, float &value, float defaultVal){
	
	bool answer = false;
	int i = 1;
	param = param + "=";
	std::string str;

	while ( (answer == false) && (i < argc) ) {
		str = argv[i];
		if (str.substr(0, param.length()) == param){
			answer = true;
			str.erase(0, param.length());
			value = std::stof(str);
		}
		i++;
	}

	if (answer == false)
		value = defaultVal;
}

void MB::readParam(int argc, char **argv, std::string param, double &value, double defaultVal){
	
	bool answer = false;
	int i = 1;
	param = param + "=";
	std::string str;

	while ( (answer == false) && (i < argc) ) {
		str = argv[i];
		if (str.substr(0, param.length()) == param){
			answer = true;
			str.erase(0, param.length());
			value = std::stod(str);
		}
		i++;
	}

	if (answer == false)
		value = defaultVal;
}

void MB::readParam(int argc, char **argv, std::string param, std::string &value, std::string defaultVal){

	bool answer = false;
	int i = 1;
	param = param + "=";
	std::string str;

	while ( (answer == false) && (i < argc) ) {
		str = argv[i];
		if (str.substr(0, param.length()) == param){
			answer = true;
			str.erase(0, param.length());
			value = str;
		}
		i++;
	}

	if (answer == false)
		value = defaultVal;
}

void MB::readParam(int argc, char **argv, std::string param, bool &value, bool defaultVal){
	
	bool answer = false;
	int i = 1;
	param = param + "=";
	std::string str;

	while ( (answer == false) && (i < argc) ) {
		str = argv[i];
		if (str.substr(0, param.length()) == param){
			answer = true;
			str.erase(0, param.length());
			if (str == "true")
				value = true;
			else
				value = false;
		}
		i++;
	}

	if (answer == false)
		value = defaultVal;
}

void MB::readParam(int argc, char **argv, std::string param, std::vector<float> &value, float defaultVal){
	
	bool answer = false;
	int i = 1;
	param = param + "=";
	std::string str;

	while ( (answer == false) && (i < argc) ) {
		str = argv[i];
		if (str.substr(0, param.length()) == param){
			answer = true;
            str.erase(0, param.length());
			std::size_t current, previous = 0;
            current = str.find(',');
            while (current != std::string::npos) {
                value.push_back(std::stof(str.substr(previous, current - previous)));
                previous = current + 1;
                current = str.find(',', previous);
            }
            value.push_back(std::stof(str.substr(previous, current - previous)));
        }
		i++;
	}

	if (answer == false)
		value.push_back(defaultVal);
}

void MB::readParam(int argc, char **argv, std::string param, std::vector<double> &value, double defaultVal){
	
	bool answer = false;
	int i = 1;
	param = param + "=";
	std::string str;

	while ( (answer == false) && (i < argc) ) {
		str = argv[i];
		if (str.substr(0, param.length()) == param){
			answer = true;
            str.erase(0, param.length());
			std::size_t current, previous = 0;
            current = str.find(',');
            while (current != std::string::npos) {
                value.push_back(std::stod(str.substr(previous, current - previous)));
                previous = current + 1;
                current = str.find(',', previous);
            }
            value.push_back(std::stod(str.substr(previous, current - previous)));
        }
		i++;
	}

	if (answer == false)
		value.push_back(defaultVal);
}