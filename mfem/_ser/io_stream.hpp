namespace PyMFEM
{
class wFILE
{
    private:
       int _isSTDOUT;
       char _filename[255];
    public:
       wFILE(void){
	 strcpy(_filename, "__stdout__");
	 _isSTDOUT = 1;					 
       }
       wFILE(const char *filename, int isSTDOUT=0){
	 strcpy(_filename, filename);
	 _isSTDOUT = isSTDOUT;
       }
       int isSTDOUT(void) {
	  return _isSTDOUT;
       }
       char * getFilename(void){
	 return _filename;
       }
};
}

