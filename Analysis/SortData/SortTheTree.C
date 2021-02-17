#include "Sort.C+"

void SortTheTree(const char *name)
{
	TChain *c = new TChain("FUEL-PIN");
	c->Add(name);

	Sort sorter(c);

	TString s1(name);
	Ssiz_t pos = s1.Last('_');
	if(pos<0) {
		cout << "must provide filename in format xxxx_t*.root" <<endl;
		return;
	}
	TString FileName = s1(0,pos+1);
	FileName.Append("sorted.root");
	sorter.Loop(FileName.Data());

}