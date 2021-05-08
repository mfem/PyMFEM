namespace PyMFEM
{
class wFILE
{
    private:
       int _isSTDOUT;
       char _filename[255];
       int _precision=8;
       bool _temporary=false;
    public:
       wFILE(void){
	 strcpy(_filename, "__stdout__");
	 _isSTDOUT = 1;
	 _temporary = 0;
       }
       wFILE(const char *filename, int precision = 8, bool temporary = 0){
	 strcpy(_filename, filename);
	 _isSTDOUT = 0;
	 _precision = precision;
	 _temporary = temporary;
       }
       int isSTDOUT(void) {
	  return _isSTDOUT;
       }
       bool isTemporary(void) {
	  return _temporary;
       }
       char * getFilename(void){
	 return _filename;
       }
       bool isGZ(void){
         int len = strlen(_filename);
         if ((len >= 3) && strcmp(_filename+len-3, ".gz") == 0){
	   return true;
	 } else {
	   return false;
	 }
       }
       int getPrecision(void){
	 return _precision;
       }
       void setPrecision(int precision){
	 _precision = precision;
       }
};
}

