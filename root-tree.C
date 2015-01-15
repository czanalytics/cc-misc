// create and analyze a ROOT tree 
// originally from http://root.cern.ch/drupal/

// Prefer compiled:
#include "EventData.h+"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

void createTree(ULong64_t numEvents = 200) {
   TFile* f = new TFile("eventdata.root", "RECREATE");
   TTree* tree = new TTree("EventTree", "Tutorial tree");
   tree->SetAutoSave(0); // to not confuse people with "EventTree;2" keys

   EventData* event = new EventData();
   tree->Branch("event", &event);

   Particle p;

   for (ULong64_t i = 0; i < numEvents; ++i) {
      event->Clear();
      int nParticles = 10 * gRandom->Exp(10.);
      if (i == 0) nParticles = 1200; // to have a large value for TSelector's array size
      for (int ip = 0; ip < nParticles; ++ip) {
         do {
            p.fPosX = gRandom->Gaus(gRandom->PoissonD(0.1), 1.)
               + gRandom->BreitWigner(0.1, 0.1);;
         } while (fabs(p.fPosX) > 10.);
         p.fPosY = gRandom->Gaus(gRandom->PoissonD(0.01), .7);
         p.fPosZ = gRandom->Gaus(gRandom->PoissonD(10), 19.);

         p.fMomentum = gRandom->Exp(12);
         p.fMomentumPhi = gRandom->Uniform(2*TMath::Pi());
         do {
            p.fMomentumEta = gRandom->BreitWigner(0.01, 10.);
         } while (fabs(p.fMomentumEta) > 12.);

         event->AddParticle(p);
      }
      event->SetSize();

      tree->Fill();

      if (i % (numEvents/50) == 0) {
         printf("*");
         fflush(stdout);
      }
   }
   printf("\n");
   tree->Write();
   delete f;
}


void AnalyzeTree()
{
   // Variables used to store the data
   TH1F *hPosX;  // X position of the particles

   // open the file
   TFile *f = TFile::Open("http://lcg-heppkg.web.cern.ch/lcg-heppkg/ROOT/eventdata.root");
   if (f == 0) {
      // if we cannot open the file, print an error message and return immediatly
      printf("Error: cannot open http://lcg-heppkg.web.cern.ch/lcg-heppkg/ROOT/eventdata.root!\n");
      return;
   }

   // Create tyhe tree reader and its data containers
   TTreeReader myReader("EventTree", f);
   // The branch "fPosX" contains doubles; access them as particlesPosX.
   TTreeReaderArray<Double_t> particlesPosX(myReader, "fParticles.fPosX");
   // The branch "fMomentum" contains doubles, too; access those as particlesMomentum.
   TTreeReaderArray<Double_t> particlesMomentum(myReader, "fParticles.fMomentum");

   // create the TH1F histogram
   hPosX = new TH1F("hPosX", "Position in X", 20, -5, 5);
   // enable bin errors:
   hPosX->Sumw2();

   // Loop over all entries of the TTree or TChain.
   while (myReader.Next()) {
      // Do the analysis...
      for (int iParticle = 0; iParticle < particlesPosX.GetSize(); ++iParticle) {
         if (particlesMomentum[iParticle] > 40.0)
            hPosX->Fill(particlesPosX[iParticle]);
      }
   }

   // Fit the histogram:
   hPosX->Fit("pol2");
   // and draw it:
   hPosX->Draw();
}
