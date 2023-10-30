# Description

GenericIO is a library tuned to seismic data, but useful in a variety of other fields,
that allows you to write IO and parameter handling using a  single interface and translates those calls in (semi-)optimal ways to different underling processing packages and hardware environments.

# Concept

For parameter handling the library supports the ability to read and write booleans, integers, floats,
doubles, and strings. Either has single valiues or comma seperated lists.

GenericIO follows the concept of SEPlib that data, even when irregular, is most easily accessed by first
describing the data as multi-dimensional hypercubes.  When performing IO the user requests sub-cubes of 
abritrary dimensions of the hypercube.  The library is responsible for fetching/putting the data to/from disk.

# Supported backends

Currently GenericIO supports completely 
 1) regular-sampled data written in SEP/Magascar format written to conventional disks
 2) regular-sampled data described by JSON files 
   a. Written as a single file to disk
   b. Written to multiple files on disk
   c. Written to GCP as objects

GenericIO partially supports trace-based data, such as standard SEG-Y and SU fomat.

Complete irregular data support is in development.

# Buffers libbrary

GenericIO uses a pacakge called buffers to handle multi-file support.  The buffers library breaks
up a multi-dimensional hypercube int sub-cubes. Each of these sub-cubes can exist in memory, compressed in memory,
or on disk (in compressed or decompressed form).  Currently only ZFP compression (up to 4-D) is supported. The
buffers library is initalized by four different objects
  1. Compresison method (How/whether to compress data)
  2. Blocking method (how to split the data into sub-cubes)
  3. Memory method (enables a caching system that will flush old sub-cubes to disk when a memory threshold has been reached)

 