//
//    AS_RENAME  rename As -> AsHypreParMatrix
//    AS_WRAP    instatiate As <T>
//
%define CONSTRUCTOR_RENAME(T)
%rename(OperatorHandleFrom ## T) mfem::OperatorHandle::OperatorHandle<mfem:: ## T >;
%enddef
%define AS_RENAME(T)
%rename(As ## T) mfem::OperatorHandle::As<mfem:: ## T >;
%enddef
%define IS_RENAME(T)
%rename(Is ## T) mfem::OperatorHandle::Is<mfem:: ## T >;
%enddef
%define GET_RENAME(T)
%rename(Get ## T) mfem::OperatorHandle::Get<mfem:: ## T >;
%enddef
%define RESET_RENAME(T)
%rename(Reset ## T) mfem::OperatorHandle::Reset<mfem:: ## T >;
%enddef
%define CONVERT_FROM_RENAME(T)
%rename(ConvertFrom ## T) mfem::OperatorHandle::ConvertFrom <mfem:: ## T >;
%enddef

%define CONSTRUCTOR_WRAP(T)
%template(OperatorHandle) mfem::OperatorHandle::OperatorHandle<T >;
%enddef
%define AS_WRAP(T)
%template(As) mfem::OperatorHandle::As<mfem:: ## T >;
%enddef
%define IS_WRAP(T)
%template(Is) mfem::OperatorHandle::Is<mfem:: ## T >;
%enddef
%define GET_WRAP(T)
%template(Get) mfem::OperatorHandle::Get<mfem:: ##T >;
%enddef
%define RESET_WRAP(T)
%template(Reset) mfem::OperatorHandle::Reset<mfem:: ## T >;
%enddef
%define CONVERT_FROM_WRAP(T)
%template(ConvertFrom) mfem::OperatorHandle::ConvertFrom<mfem:: ## T >;
%enddef
