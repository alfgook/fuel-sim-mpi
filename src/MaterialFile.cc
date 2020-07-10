#include "MaterialFile.hh"

#include <algorithm>
#include <fstream>
#include <functional>

#include "G4Material.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4SystemOfUnits.hh"

struct MaterialFileEntry {
	G4int A;
	G4int Z;
	G4String element;
	G4String nuclide;
	G4double atomic_mass;
	G4double mass;
	G4double moles;
	G4int ZAID;
	G4double abundance;

	bool operator > (const MaterialFileEntry& str) const
    {
        if(Z == str.Z) return (A > str.A);
        return (Z > str.Z);
    }
	bool operator < (const MaterialFileEntry& str) const
    {
        if(Z == str.Z) return (A < str.A);
        return (Z < str.Z);
    }
    bool operator == (const MaterialFileEntry& str) const
    {
    	return (Z == str.Z);
    }
};

MaterialFile::MaterialFile(const char *FileName, G4String matName, G4double Volume)
{
	fVolume = Volume;

	// a list of important (alpha,n) target nuclei according to Table 5.2.10 of Scale 6.2 manual
	std::array<G4int,19> alphaNtargets = {30070,40090,50100,50110,60130,70140,80170,80180,90190,100210,100220,110230,120250,120260,130270,140290,140300,150310,170370};
	G4cout << "======================================================" << G4endl;
	G4cout << "== Reading material from : " << FileName << G4endl;
	G4cout << "======================================================" << G4endl;

		//read AME-2016
	std::string str;
	std::ifstream textfile("input/atomic-mass-eval-2016.txt");
	if(!textfile.is_open()) {
	  G4ExceptionDescription ed;
	  ed << "      Could not open the file: input/input/atomic-mass-eval-2016.txt" << G4endl;
	  G4Exception("MaterialFile::MaterialFile(...)",
	              "MaterialFile.01",
	              FatalException,
	              ed);
	}
	getline(textfile,str);
	getline(textfile,str);

	std::vector<int> AME_Z;
	std::vector<int> AME_A;
	std::vector<double> AME_mass;
	while(!textfile.eof()) {
		int Z, N, A;
		double mass;
		textfile >> N >> Z >> A >> mass;
		AME_Z.push_back(Z);
		AME_A.push_back(A);
		AME_mass.push_back(mass*1.e-06);
	}
	textfile.close();

	//--read file used to convert element name to Z
	std::vector<std::string> elements;

	textfile.open("input/elements.txt");
	if(!textfile.is_open()) {
		G4ExceptionDescription ed;
		ed << "      Could not open the input file: input/elements.txt" << G4endl;
		G4Exception("MaterialFile::MaterialFile(...)",
		            "MaterialFile.02",
		            FatalException,
		            ed);
	}
	while(!textfile.eof()) {
		getline(textfile,str);
		elements.push_back(str);
	}
	textfile.close();

	//read the plt concentration file

	textfile.open(FileName);
	if(!textfile.is_open()) {
		G4ExceptionDescription ed;
		ed << "      Could not open the input file: " << FileName << G4endl;
		G4Exception("MaterialFile::MaterialFile(...)",
		            "MaterialFile.03",
		            FatalException,
		            ed);
	}


	for(int i=0;i<6;i++) getline(textfile,str); //skip header
	std::streampos pos = textfile.tellg();

	//count the lines
	G4int number_of_lines = 0;
    while (std::getline(textfile, str)) {
    	if(str.find("subtotal")!=std::string::npos) break;
    	++number_of_lines;
    }
	theList.reserve(number_of_lines-2);

	textfile.seekg(pos); //rewind the stream

	std::string nuclide;
	G4double mass;
	total_mass = 0;
	while(!textfile.eof()) {
		textfile >> nuclide >> mass;
		if(nuclide=="subtotal") break;
		mass *= g;
		total_mass += mass;

		MaterialFileEntry anEntry;
		anEntry.mass = mass;
		anEntry.nuclide = nuclide;

		std::string element;
    	element.reserve( nuclide.size() );
    	std::string Astring;
    	Astring.reserve( 3 );

    	for ( char c : nuclide ) {
	        if ( std::isdigit( c ) ) {
	        	Astring.push_back(c);
	        }
	    }

    	for ( char c : nuclide ) {
	        if ( std::isdigit( c ) ) {
	        	break;
	        }
	        element.push_back(c);
	    }

	    anEntry.element = element;
	    anEntry.A = atoi(Astring.data());
		for(size_t i=0;i<elements.size();i++) {
			if(element==elements[i]) {
				anEntry.Z = i + 1;
				break;
			}
		}

		for(size_t k=0;k<AME_Z.size();k++) {
			if(anEntry.Z==AME_Z[k] && anEntry.A==AME_A[k]) {
				anEntry.atomic_mass = AME_mass[k]*g/mole;
			}
		}

		anEntry.moles = mass/anEntry.atomic_mass;
		anEntry.ZAID = anEntry.Z*10000 + anEntry.A*10;
		//anEntry.ZAID = anEntry.Z*1000 + anEntry.A;

		//G4cout << nuclide << "  " << anEntry.A << " " << anEntry.Z << " " << anEntry.element << G4endl;
		theList.push_back(anEntry);

	}
	textfile.close();
	//std::sort(theList.begin(),theList.end(),std::greater<MaterialFileEntry>());
	std::sort(theList.begin(),theList.end());

	//add up and remove duplicates
	for(size_t i=0;i<theList.size();i++) {
		if(theList[i]==theList[i+1]) {
			if(theList[i].A==theList[i+1].A) {
				theList[i].mass += theList[i+1].mass;
				theList.erase(theList.begin()+(i+1));
			}
		}
	}

	//remove entries with small mass
	G4int nn = 0;
	for(size_t i=0;i<theList.size();i+=nn) {

		nn = 1;
		//never remove the (alpha,n) targets
		if( std::find(std::begin(alphaNtargets),std::end(alphaNtargets),theList[i].ZAID) != std::end(alphaNtargets) ) {
			continue;
		}

		if(theList[i].mass<0.1*g) {
			//G4cout << "removing " << theList[i].nuclide << " with mass " << theList[i].mass/g << ", nn = " << nn << G4endl;
			//G4cout << "------" << G4endl;
			theList.erase(theList.begin()+i,theList.begin()+(i+nn));
			nn = 0;
		}
	}

	for(size_t i=0;i<theList.size();i++) {
		G4int Z = theList[i].Z;

		G4double sum_of_moles = 0.;
		G4int nComp = 0;
		for(size_t j=i;j<theList.size();j++) {
			if(theList[j].Z!=Z) break; 
			sum_of_moles += theList[j].moles;
			++nComp;
		}

		for(size_t j=i;j<theList.size();j++) {
			if(theList[j].Z!=Z) break;
			theList[j].abundance = theList[j].moles/sum_of_moles;
			if(theList[j].abundance<1.E-04) theList[j].Z = 0; //set Z=0 as a flag to remove the line
		}
	}

	Print();

	MaterialFileEntry temp;
	temp.Z = 0;
	theList.erase(std::remove(std::begin(theList),std::end(theList),temp),theList.end());
	Print();
	//Print();
/*
	// remove elements with small mass fraction
	G4int nn = 0;
	for(G4int i=0;i<theList.size();i+=nn) {

		nn = 0;
		G4double mass_of_element = 0;
		for(G4int j=i;j<theList.size();j++) {
			if(theList[j].Z!=theList[i].Z) break;
			++nn;
			mass_of_element += theList[j].mass;
		}

		//if(mass_of_element/total_mass<1.E-09) {
		if(mass_of_element<1.*g) {
			G4cout << "removing " << theList[i].nuclide << " with mass " << mass_of_element/g << ", nn = " << nn << G4endl;
			G4cout << "removing " << theList[i+nn].nuclide << " with mass " << mass_of_element/g << G4endl;
			G4cout << "------" << G4endl;
			theList.erase(theList.begin()+i,theList.begin()+(i+nn));
			nn = 0;
		}
	}

	total_mass = 0;
	for(G4int i=0;i<theList.size();i++) {
		total_mass += theList[i].mass;
	}
	Print();
	G4cout << "removed elements with small mass" << G4endl;


	// remove isotopes with small abundance
	for(G4int i=0;i<theList.size();i++) {
		G4double sum_of_moles = 0.;
		for(G4int j=i;j<theList.size();j++) {
			if(theList[j].Z!=theList[i].Z) break; 
			sum_of_moles += theList[j].moles;
		}

		for(G4int j=i;j<theList.size();j++) {
			if(theList[j].Z!=theList[i].Z) break; 
			G4double abundance = theList[j].moles/sum_of_moles;
			if(abundance<1.e-4) {
				theList.erase(theList.begin()+j,theList.begin()+(j+1));
			}
		}
	}
	Print();
	G4cout << "removed isotopes with small abundance" << G4endl;
*/
	total_mass = 0.;
	nElements = 1;
	G4int ZZ = theList[0].Z;
	for(size_t i=0;i<theList.size();i++) {
		if(theList[i].Z!=ZZ) {
			nElements++;
			ZZ = theList[i].Z;
		}
		total_mass += theList[i].mass;
	}
	PrintZAID();

	/*G4cout << "nElements = " << nElements << G4endl;
	G4cout << "matName = " << matName << G4endl;*/

	G4double density = total_mass/Volume;
	theMaterial = new G4Material(matName,density,nElements,kStateSolid);
	G4int nComp = 0;
	for(size_t i=0;i<theList.size();i+=nComp) {
		G4cout << "------------------"  << G4endl;
		G4int Z = theList[i].Z;

		G4double sum_of_moles = 0.;
		nComp = 0;
		for(size_t j=i;j<theList.size();j++) {
			if(theList[j].Z!=Z) break; 
			sum_of_moles += theList[j].moles;
			++nComp;
		}

		G4double mass_of_element = 0;
		for(size_t j=i;j<theList.size();j++) {
			if(theList[j].Z!=Z) break;
			mass_of_element += theList[j].mass;
		}

		//G4cout << "Element = " <<  theList[i].element << G4endl;
		//G4cout << "nComp = " <<  nComp << G4endl;
		G4Element *anElement = new G4Element(theList[i].element,theList[i].element,nComp);

		for(size_t j=i;j<theList.size();j++) {
			if(theList[j].Z!=Z) break;
			G4double abundance = theList[j].moles/sum_of_moles;
			G4double A = theList[j].A;
			G4double atomic_mass = theList[j].atomic_mass;
			anElement->AddIsotope(new G4Isotope(theList[j].nuclide,Z,A,atomic_mass),abundance);

			/*G4cout << "----"  << G4endl;
			G4cout << "  isotope = " << theList[j].nuclide << G4endl;
			G4cout << "  Z = " << Z << G4endl;
			G4cout << "  A = " << A << G4endl;
			G4cout << "  atomic_mass = " << atomic_mass/(g/mole) << G4endl;
			G4cout << "  abundance = " << abundance << G4endl;*/
		}

		theMaterial->AddElement(anElement,mass_of_element/total_mass);
	}
}

MaterialFile::~MaterialFile()
{

}

void MaterialFile::Print()
{
	G4cout << "------MaterialFile::Print()------------------" << G4endl;
	for(size_t i=0;i<theList.size();i++) {
		G4cout << theList[i].nuclide << " " << theList[i].Z << " " << theList[i].A << " " << theList[i].element << " " << theList[i].mass/g << G4endl;
	}
	G4cout << "tot mass = " << total_mass/g << G4endl;
	G4cout << "density = " << (total_mass/g)/(fVolume/cm3) << G4endl;
	G4cout << "entries = " << theList.size() << G4endl;
	G4cout << "---------------------------------------------" << G4endl;
}

void MaterialFile::PrintZAID()
{
	G4cout << "------MaterialFile::PrintZAID()------------------" << G4endl;
	for(size_t i=0;i<theList.size();i++) {
		//G4cout << theList[i].nuclide << " " << theList[i].ZAID << " " << theList[i].Z << " " << theList[i].A << " " << theList[i].mass/total_mass << G4endl;
		G4int tmp = theList[i].Z*1000 + theList[i].A;
		G4cout.precision(9);
		G4cout << tmp << ".70c -" << std::scientific << theList[i].mass/total_mass << G4endl;
	}
	G4cout << "tot mass = " << total_mass/g << G4endl;
	//G4cout << "density = " << (total_mass/g)/(fVolume/cm3) << G4endl;
	G4cout << "entries = " << theList.size() << G4endl;
	G4cout << "---------------------------------------------" << G4endl;
}

G4Material* MaterialFile::GetMaterial()
{
	return theMaterial;
}
