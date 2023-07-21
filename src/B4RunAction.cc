//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B4RunAction.cc
/// \brief Implementation of the B4RunAction class

#include "B4RunAction.hh"
#include "B4Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4RunAction::B4RunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(0);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in B4Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  ///G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(0);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple
  //
  
  // Creating histograms
  analysisManager->CreateH1("Ecore","Edep in core", 100, 0., 800*MeV);
  analysisManager->CreateH1("Eann1","Edep in ann1", 100, 0., 100*MeV);
  analysisManager->CreateH1("Lcore","trackL in core", 100, 0., 1*m);
  analysisManager->CreateH1("Lann1","trackL in ann1", 100, 0., 50*cm);

  analysisManager->CreateH1("Eann2","Edep in ann2", 100, 0., 800*MeV);
  analysisManager->CreateH1("Eann3","Edep in ann3", 100, 0., 800*MeV);
  analysisManager->CreateH1("Eann4","Edep in ann4", 100, 0., 800*MeV);
  analysisManager->CreateH1("Eann5","Edep in ann5", 100, 0., 800*MeV);
  analysisManager->CreateH1("Eann6","Edep in ann6", 100, 0., 800*MeV);
  // Creating ntuple
  //
  analysisManager->CreateNtuple("B4", "Edep and TrackL");
  analysisManager->CreateNtupleDColumn("Ecore");
  analysisManager->CreateNtupleDColumn("Eann1");
  analysisManager->CreateNtupleDColumn("Lcore");
  analysisManager->CreateNtupleDColumn("Lann1");
  
  analysisManager->CreateNtupleDColumn("Eann2");
  analysisManager->CreateNtupleDColumn("Eann3");
  analysisManager->CreateNtupleDColumn("Eann4");
  analysisManager->CreateNtupleDColumn("Eann5");
  analysisManager->CreateNtupleDColumn("Eann6");
  
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4RunAction::~B4RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "RUN_RESULTS";
  analysisManager->OpenFile(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->GetH1(1) ) {
    G4cout << G4endl << " ----> print histograms statistic ";
    if(isMaster) {
      G4cout << "for the entire run " << G4endl << G4endl; 
    }
    else {
      G4cout << "for the local thread " << G4endl << G4endl; 
    }
    
    G4cout << " Ecore : mean = " 
       << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy") 
       << " rms = " 
       << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;
    
    G4cout << " Eann1 : mean = " 
       << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy") 
       << " rms = " 
       << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;
    
    G4cout << " Lcore : mean = " 
      << G4BestUnit(analysisManager->GetH1(2)->mean(), "Length") 
      << " rms = " 
      << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Length") << G4endl;

    G4cout << " Lann1 : mean = " 
      << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length") 
      << " rms = " 
      << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Length") << G4endl;

    G4cout << " Eann2 : mean = "
       << G4BestUnit(analysisManager->GetH1(4)->mean(), "Energy")
       << " rms = "
       << G4BestUnit(analysisManager->GetH1(4)->rms(),  "Energy") << G4endl;

    G4cout << " Eann3 : mean = "
       << G4BestUnit(analysisManager->GetH1(5)->mean(), "Energy")
       << " rms = "
       << G4BestUnit(analysisManager->GetH1(5)->rms(),  "Energy") << G4endl;

    G4cout << " Eann4 : mean = "
       << G4BestUnit(analysisManager->GetH1(6)->mean(), "Energy")
       << " rms = "
       << G4BestUnit(analysisManager->GetH1(6)->rms(),  "Energy") << G4endl;
	        
    G4cout << " Eann5 : mean = "
       << G4BestUnit(analysisManager->GetH1(7)->mean(), "Energy")
       << " rms = "
       << G4BestUnit(analysisManager->GetH1(7)->rms(),  "Energy") << G4endl;

    G4cout << " Eann6 : mean = "
       << G4BestUnit(analysisManager->GetH1(8)->mean(), "Energy")
       << " rms = "
       << G4BestUnit(analysisManager->GetH1(8)->rms(),  "Energy") << G4endl;
 }
  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
