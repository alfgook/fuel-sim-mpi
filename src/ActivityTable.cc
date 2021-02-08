#include "Config.h"
#include "ActivityTable.hh"
#include "Randomize.hh"
#include "G4DecayTable.hh"
#include "G4IonTable.hh"

/*#ifdef G4MULTITHREADED
#include "G4AutoLock.hh"
G4Mutex ActivityTable::ActivityTableMutex = G4MUTEX_INITIALIZER; //used for reading from the input text file
#endif*/

#include "G4MPImanager.hh"

ActivityTable::ActivityTable(G4String aFile, MyRadioactiveDecayBase *aRadDecay)
{
  fBin = 0;
  fRadDecay = aRadDecay;
  /*#ifdef G4MULTITHREADED
  G4AutoLock lock(&ActivityTable::ActivityTableMutex); //lock the mutex while reading from the text-file
  #endif*/

  //fRadDecay->LoadAllDecayTables();

  G4String dir(INPUT_DIR);
  if(!aFile.size()) aFile = dir + "/activities.txt";
  fInit = true;

  nPrimaries = 0;
  //---read the file-----------------------------------
  std::string str;
  std::ifstream textfile(aFile.data());
  if(!textfile.is_open()) {
    G4ExceptionDescription ed;
      ed << "      Could not open activity file: " << aFile.data() << G4endl;
      G4Exception("ActivityTable::ActivityTable(...)",
                  "ActivityTable.01",
                  FatalException,
                  ed);

    /*textfile.open("input/DummyActivity.txt");
    if(!textfile.is_open()) {
      G4ExceptionDescription ed2;
      ed2 << "      Could not open the dummy activity file: input/DummyActivity.txt" << G4endl;
      G4Exception("ActivityTable::ActivityTable(...)",
                  "ActivityTable.02",
                  FatalException,
                  ed2);
    }*/

    fInit = false;
  }

  for(G4int i=0;i<7;i++) getline(textfile,str);

  std::vector<std::string> nuclides;
  std::vector<double> activity1;
  std::vector<int> BRbias1;
  while(!textfile.eof()) {
    getline(textfile,str);
    if(str.find("------------")!=std::string::npos) break;
    std::stringstream ss(str);
    std::string s1;
    double act;
    int bb = 0;
    ss >> s1 >> act >> bb;
    nuclides.push_back(s1);
    activity1.push_back(act);
    BRbias1.push_back(bb);
  }
  textfile.close();
  //---end read the file-----------------------------------

  //--read file used to convert element name to Z
  std::vector<std::string> elements;

  aFile = dir + "/elements.txt";
  textfile.open(aFile.data());
  if(!textfile.is_open()) {
    G4ExceptionDescription ed;
    ed << "      Could not open the input file: input/elements.txt" << G4endl;
    G4Exception("ActivityTable::ActivityTable(...)",
                "ActivityTable.03",
                FatalException,
                ed);
  }
  while(!textfile.eof()) {
    getline(textfile,str);
    elements.push_back(str);
  }
  textfile.close();

  std::vector<std::string> ExcitationKey;
  std::vector<double> ExcitationEnergy;

  aFile = dir + "/MetaStables.txt";
  textfile.open(aFile.data());
  if(!textfile.is_open()) {
    G4ExceptionDescription ed;
    ed << "      Could not open the input file: input/MetaStables.txt" << G4endl;
    G4Exception("ActivityTable::ActivityTable(...)",
                "ActivityTable.04",
                FatalException,
                ed);
  }
  getline(textfile,str);
  while(!textfile.eof()) {
    getline(textfile,str);
    if(str[0]=='#') continue;
    std::stringstream ss(str);
    std::string s1;
    double ee;
    ss >> s1 >> ee;
    ExcitationKey.push_back(s1);
    ExcitationEnergy.push_back(ee);
  }
  textfile.close();

  fZZ.clear();
  fAA.clear();
  activity.clear();
  fBRbias.clear();
  fTables.clear();
  for(unsigned int i=0;i<nuclides.size();i++) {

    int A = -1;
    int Z = -1;
    bool BRbias = false;
    bool MetaStable = false;
    G4double ExcEnergy = 0.;

    //Including the meta-stable states
    size_t pos = nuclides[i].find("-");
    std::string element = nuclides[i].substr(0,pos);
    for(unsigned int j=0;j<elements.size();j++) {
      if(elements.at(j)==element) {
        Z = j+1;
        break;
      }
    }
    size_t pos2 = nuclides[i].size();
    if(nuclides[i].back()=='m') {
      --pos2;
      MetaStable = true;
      G4bool found = false;
      for(unsigned int k=0;k<elements.size();k++) {
        if(ExcitationKey[k]==nuclides[i]) {
          ExcEnergy = ExcitationEnergy[k]/1000.;
          found = true;
          break;
        }
      }
      if(!found) {
        ExcEnergy = 0.;
        G4ExceptionDescription ed;
        ed << "      could not find an entry for " << nuclides[i] << " in the file input/MetaStables.txt " <<
              "the excitation energy will be set to 0!!! You can add an entry into the file to get rid of this message" << G4endl;
        G4Exception("ActivityTable::ActivityTable(...)",
                  "ActivityTable.05",
                  JustWarning,
                  ed);
      }
    } else {
      MetaStable = false;
      ExcEnergy = 0.;
    }
    std::string strA = nuclides[i].substr(pos+1,pos2-(pos+1));
    A = atoi(strA.data());

    if(BRbias1[i]) {
      BRbias = true;
    } else {
      BRbias = false;
    }

    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,ExcEnergy);
    if(!ion) continue;
    G4DecayTable *theDecayTable = fRadDecay->GetDecayTable(ion);
    if(!theDecayTable) continue;

    activity.push_back(activity1[i]);
    fAA.push_back(A);
    fZZ.push_back(Z);
    fBRbias.push_back(BRbias);
    fMetaStable.push_back(MetaStable);
    fExcEnergy.push_back(ExcEnergy);
    //fTables.push_back(MyDecayTable table);
//    G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << nuclides[i] << "  " << Z << "  " << A << "  " << fExcEnergy[i] << G4endl;
  }
  //--------------------------------------------------------

  //print the activity table to screen
  /*for(unsigned int i=0;i<activity.size();i++) {
    G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << fZZ[i] << "\t" << fAA[i] << "\t" << activity[i] << "\t" << fExcEnergy[i] << G4endl;
  }*/

  //std::vector<double> activityCumulative;
  activityCumulative.clear();
  activityCumulative.push_back(activity[0]);
  for(unsigned int i=1;i<activity.size();i++) activityCumulative.push_back(activity[i]+activityCumulative[i-1]);

  activityTotal = activityCumulative.back();
  for(unsigned int i=0;i<activity.size();i++) activityCumulative[i] /= activityTotal;

//  G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "Activity = " << activityTotal << G4endl; 
}

ActivityTable::~ActivityTable()
{

}

void ActivityTable::RestrictTo(G4String KinematicsName)
{
    /*#ifdef G4MULTITHREADED
    G4AutoLock lock(&ActivityTable::ActivityTableMutex); //lock the mutex while reading from the text-file
    #endif*/
    G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "====================================" << G4endl;
    G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "==  Restricting decay to " << KinematicsName << G4endl;
    G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "==  total activity = " << activityTotal << G4endl;

    G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "   Z   A   rest-activity  total-activity" << G4endl;

  for(size_t bin=0;bin<activityCumulative.size();bin++) {
    G4int Z = fZZ.at(bin);
    G4int A = fAA.at(bin);
    G4double E = fExcEnergy.at(bin);

    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,E);
    if(!ion) continue;
    G4DecayTable *theDecayTable = fRadDecay->GetDecayTable(ion);
    if(!theDecayTable) continue;

    G4int ndecaych = 0;
    G4double BRsum = 0.;
    MyDecayTable table;
    for(G4int i=0;i<theDecayTable->entries();i++) {
      if(theDecayTable->GetDecayChannel(i)->GetKinematicsName()==KinematicsName) {
        ndecaych++;
        BRsum += theDecayTable->GetDecayChannel(i)->GetBR();
      }
    }
    table.SetBrsum(BRsum);
    for(G4int i=0;i<theDecayTable->entries();i++) {
      if(theDecayTable->GetDecayChannel(i)->GetKinematicsName()==KinematicsName) {
        G4double BR = theDecayTable->GetDecayChannel(i)->GetBR()/BRsum;
        table.Add(BR,i);
      }
    }
    if(!ndecaych) { //no SF decay channel for this nuclide
      //G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "no " << KinematicsName << " for (Z,A) = (" << fZZ[bin] <<"," << fAA[bin] << ")" << G4endl;
      activity[bin] = 0;
    } else {
      G4double BR = table.GetBrsum();
      G4double newActivity = BR*activity.at(bin);
      activity[bin] = newActivity;
      G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "   " << fZZ[bin] << "  " << fAA[bin] << "  " << newActivity <<  "  " << newActivity/BR << G4endl;
      for(size_t i=0;i<table.GetEntries();i++) G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "   " << table.GetEntry(i);
      G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << G4endl;
    }
    fTables.push_back(table);
  }

  activityCumulative.clear();
  activityCumulative.push_back(activity[0]);
  for(unsigned int i=1;i<activity.size();i++) activityCumulative.push_back(activity[i]+activityCumulative[i-1]);

  activityTotal = activityCumulative.back();
  for(unsigned int i=0;i<activity.size();i++) activityCumulative[i] /= activityTotal;

  G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "== restricted activity = " << activityTotal << G4endl; 
  G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "====================================" << G4endl;

  CleanUpTable();
  /*#ifdef G4MULTITHREADED
  lock.unlock(); //explicit unlock
  #endif*/
}

void ActivityTable::ExcludeAphaAndSF()
{
    /*#ifdef G4MULTITHREADED
    G4AutoLock lock(&ActivityTable::ActivityTableMutex); //lock the mutex while reading from the text-file
    #endif*/
    G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "====================================" << G4endl;
    G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "==  Restricting decay to exclude (sf) and alpha decay" << G4endl;
    G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "==  total activity = " << activityTotal << G4endl;

    G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "   Z   A   rest-activity  total-activity  main-branch" << G4endl;

  for(size_t bin=0;bin<activityCumulative.size();bin++) {
    G4int Z = fZZ.at(bin);
    G4int A = fAA.at(bin);
    G4double E = fExcEnergy.at(bin);

    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,E);
    if(!ion) {
      G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "!ion" << G4endl;
      continue;
    }
    G4DecayTable *theDecayTable = fRadDecay->GetDecayTable(ion);
    if(!theDecayTable) {
      G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "!ion" << G4endl;
      continue;
    }
    //if(theDecayTable->entries()) theDecayTable->DumpInfo();

    G4int ndecaych = 0;
    G4double BRsum = 0.;
    MyDecayTable table;
    G4String MainBranchName;
    G4double BRmax = 0;
    for(G4int i=0;i<theDecayTable->entries();i++) {
      G4String kinName = theDecayTable->GetDecayChannel(i)->GetKinematicsName();
      if(kinName != "alpha decay" && kinName != "SF decay") {
        ndecaych++;
        BRsum += theDecayTable->GetDecayChannel(i)->GetBR();
        if(theDecayTable->GetDecayChannel(i)->GetBR()>BRmax) {
          BRmax = theDecayTable->GetDecayChannel(i)->GetBR();
          MainBranchName = kinName;
        }
      }
    }
    table.SetBrsum(BRsum);
    for(G4int i=0;i<theDecayTable->entries();i++) {
      G4String kinName = theDecayTable->GetDecayChannel(i)->GetKinematicsName();
      if(kinName != "alpha decay" && kinName != "SF decay") {
        G4double BR = theDecayTable->GetDecayChannel(i)->GetBR()/BRsum;
        table.Add(BR,i);
      }
    }
    if(!ndecaych) { //only alpa- and SF- decay channels for this nuclide
      //G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "no " << KinematicsName << " for (Z,A) = (" << fZZ[bin] <<"," << fAA[bin] << ")" << G4endl;
      activity[bin] = 0;
    } else {
      G4double BR = table.GetBrsum();
      G4double newActivity = BR*activity.at(bin);
      activity[bin] = newActivity;
      G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "   "<< fZZ[bin] << "  " << fAA[bin] << "  " << newActivity <<  "  " << newActivity/BR << "  " << MainBranchName << G4endl;
      //for(size_t i=0;i<table.GetEntries();i++) G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "   " << table.GetEntry(i);
      //G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << G4endl;
    }
    fTables.push_back(table);
  }

  CleanUpTable();

  activityCumulative.clear();
  activityCumulative.push_back(activity[0]);
  for(unsigned int i=1;i<activity.size();i++) activityCumulative.push_back(activity[i]+activityCumulative[i-1]);

  activityTotal = activityCumulative.back();
  for(unsigned int i=0;i<activity.size();i++) activityCumulative[i] /= activityTotal;

  G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "== restricted activity = " << activityTotal << G4endl; 
  G4cout<<"MPIrank"<<G4MPImanager::GetManager()->GetRank()<<" : " << "====================================" << G4endl;

  /*#ifdef G4MULTITHREADED
  lock.unlock(); //explicit unlock
  #endif*/
}

void ActivityTable::CleanUpTable()
{
  G4double max = 0.;
  for(auto it=activity.begin();it!=activity.end();it++) {
    if(max<*it) max = *it;
  }
  for(size_t bin=0;bin<activity.size();bin++) {
    G4double relative = activity.at(bin)/max;
    //if(activity.at(bin)==0) {
    if(relative<1.E-15) {
      activity.erase(activity.begin() + bin);
      activityCumulative.erase(activityCumulative.begin() + bin);
      fZZ.erase(fZZ.begin() + bin);
      fAA.erase(fAA.begin() + bin);
      fBRbias.erase(fBRbias.begin() + bin);
      fMetaStable.erase(fMetaStable.begin() + bin);
      fExcEnergy.erase(fExcEnergy.begin() + bin);
      fTables.erase(fTables.begin() + bin);
      --bin;
    }
  }
}

void ActivityTable::GenerateMotherNuclide(G4int &Z, G4int &A, G4double &E, G4bool &ApplyBRbias)
{
  G4double urn = G4UniformRand();
  for(fBin=0;fBin<activityCumulative.size();fBin++) {
    if(urn<=activityCumulative.at(fBin)) break;
  }

  try { Z = fZZ.at(fBin); }
  catch (const std::out_of_range& oor) {
    G4cerr << "GenerateMotherNuclide fZZ.at(fBin)" << G4endl;
  }
  try {A = fAA.at(fBin);}
  catch (const std::out_of_range& oor) {
    G4cerr << "GenerateMotherNuclide fAA.at(fBin)" << G4endl;
  }
  try {E = fExcEnergy.at(fBin);}
  catch (const std::out_of_range& oor) {
    G4cerr << "GenerateMotherNuclide fExcEnergy.at(fBin)" << G4endl;
  }
  try {ApplyBRbias = fBRbias.at(fBin);}
  catch (const std::out_of_range& oor) {
    G4cerr << "GenerateMotherNuclide fBRbias.at(fBin)" << G4endl;
  }

}