#include <stdio.h>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
namespace io = boost::iostreams;
typedef io::stream_buffer<io::file_descriptor_sink> boost_ofdstream;
typedef io::stream_buffer<io::file_descriptor_source> boost_ifdstream;
