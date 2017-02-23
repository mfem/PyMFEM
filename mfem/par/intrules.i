%module intrules
%{
#include "fem/intrules.hpp"
%}
%import "array.i"

%immutable IntRules;
%immutable RefinedIntRules;


%ignore mfem::Array<mfem::IntegrationPoint>::Union;
%ignore mfem::Array<mfem::IntegrationPoint>::Find;
%ignore mfem::Array<mfem::IntegrationPoint>::Sort;
%ignore mfem::Array<mfem::IntegrationPoint>::DeleteFirst;
%ignore mfem::Array<mfem::IntegrationPoint>::Unique;
%ignore mfem::Array<mfem::IntegrationPoint>::PartialSum;
%ignore mfem::Array<mfem::IntegrationPoint>::Sum;
%ignore mfem::Array<mfem::IntegrationPoint>::IsSorted;
%ignore mfem::Array<mfem::IntegrationPoint>::Save;
%ignore mfem::Array<mfem::IntegrationPoint>::Max;
%ignore mfem::Array<mfem::IntegrationPoint>::Min;
%ignore mfem::Array<mfem::IntegrationPoint>::Print;
%ignore mfem::Array<mfem::IntegrationPoint>::Load;
%template(IntegrationPointArray) mfem::Array<mfem::IntegrationPoint>;

%include "fem/intrules.hpp"



 // following symbos are ignored to wrap Array<IntegrationPoint>
 /*
%ignore Unique;
%ignore Sort;
%ignore Union;
%ignore Find;
%ignore DeleteFirst;
%ignore PartialSum;
%ignore Sum;
%ignore Save;
%ignore Print;
%ignore IsSorted;
//%namespace mfem{
%template(IntegrationPointArray) mfem::Array<mfem::IntegrationPoint>;
//}

*/
