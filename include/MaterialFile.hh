

#ifndef MaterialFile_h
#define MaterialFile_h 1

#include <vector>
#include "globals.hh"

struct MaterialFileEntry;
class G4Material;

class MaterialFile {
  public:
  	MaterialFile(const char*, G4String, G4double);
  	virtual ~MaterialFile();

  	void Print();
  	void PrintZAID();
  	G4Material* GetMaterial();

  private:
  	std::vector<MaterialFileEntry> theList;
  	G4Material *theMaterial;
  	G4double total_mass;
  	G4double fVolume;
  	G4int nElements;
};

#endif