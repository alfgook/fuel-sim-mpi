
#ifndef ActivityTable_h
#define ActivityTable_h 1

#include "MyRadioactiveDecayBase.hh"

#include "globals.hh"
#include <vector>
#include <fstream>
#include "Randomize.hh"

class MyRadioactiveDecayBase;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class MyDecayTable
{
  public:
    MyDecayTable() {
        nChannels = 0;
    };
    ~MyDecayTable() {};

    void Add(G4double aBR,G4int aNbr) {
        BR.push_back(aBR);
        ChannelNbr.push_back(aNbr);
        ++nChannels;
    }
    void SetBrsum(G4double aBRsum) { BRsum = aBRsum; }
    G4double GetBrsum() { return BRsum; }
    G4int SelectDecayChannel() {
        if(nChannels==1) return ChannelNbr[0];
        G4double rand = G4UniformRand();
        G4double sum = 0.;
        for(G4int i=0;i<nChannels;i++) {
            sum += BR[i];
            if(rand<sum) return ChannelNbr[i];
        }
        return ChannelNbr[nChannels-1];
    }
    size_t GetEntries() const { return ChannelNbr.size(); }
    G4int GetEntry(size_t i) const {
        if(i<ChannelNbr.size()) return ChannelNbr.at(i);
        return -1;
    }

  private:
    G4int                 nChannels;
    std::vector<G4double> BR;
    std::vector<G4int>    ChannelNbr;
    G4double              BRsum;

};

class ActivityTable
{
  public:
    ActivityTable(G4String,MyRadioactiveDecayBase*);    
   ~ActivityTable();

  public:
    void GenerateMotherNuclide(G4int &Z, G4int &A, G4double &E, G4bool &ApplyBRbias);
    G4double GetTotalActivity() { return activityTotal; };
    G4bool IsInit() { return fInit; };
    G4bool GetMetaStable() { 
        try {fMetaStable.at(fBin)}
        catch (const std::out_of_range& oor) {
            G4cerr << "GetMetaStable fMetaStable.at(fBin)" << G4endl;
        }
        return fMetaStable.at(fBin); 
    };
    G4int GetRestrictedDecayChannelNbr() {
        try {fTables.at(fBin)}
        catch (const std::out_of_range& oor) {
            G4cerr << "GetRestrictedDecayChannelNbr fTables.at(fBin)" << G4endl;
        }
        return fTables.at(fBin).SelectDecayChannel(); 

    };

    void RestrictTo(G4String);
    void ExcludeAphaAndSF();
    void CleanUpTable();
            
  private:

    unsigned int fBin;
    G4bool fInit;
    G4int nPrimaries;
    double activityTotal;
    std::vector<double> activity;
    std::vector<double> activityCumulative;
    std::vector<int> fZZ;
    std::vector<int> fAA;
    std::vector<bool> fBRbias;
    std::vector<bool> fMetaStable;
    std::vector<double> fExcEnergy;
    //keeps track of the channel number if the decay have been restricted to a particular channel
    std::vector<MyDecayTable> fTables;

    MyRadioactiveDecayBase *fRadDecay;

    #ifdef G4MULTITHREADED
  public:
    static G4Mutex ActivityTableMutex;
  #endif

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
