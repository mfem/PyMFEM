namespace PyMFEM
{
class wFILE
{
    private:
       int _isSTDOUT;
       char _filename[255];
       int _precision=8;
    public:
       wFILE(void){
	 strcpy(_filename, "__stdout__");
	 _isSTDOUT = 1;					 
       }
       wFILE(const char *filename, int precision = 8){
	 strcpy(_filename, filename);
	 _isSTDOUT = 0;
	 _precision = precision;
       }
       int isSTDOUT(void) {
	  return _isSTDOUT;
       }
       char * getFilename(void){
	 return _filename;
       }
       int getPrecision(void){
	 return _precision;
       }
       void setPrecision(int precision){
	 _precision = precision;
       }
};
}

