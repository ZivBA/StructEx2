#include "Vector3.h"
#include "Atom.h"
#include "RigidTrans3.h"
#include "Matrix3.h"
#include "Molecule.h"
#include "PDB.h"
#include "Match.h"
#include "GeomHash.h"
#include "Timer.h"

#include <iostream>

int oldMain(int argc , char* argv[]){
  Timer overallTime("Overall Run Time");
  overallTime.start();

  if(argc !=5) {
    cerr << "Usage: "<<argv[0]<< " target.pdb model.pdb num_rotations dist_threshold" << endl;
    exit(1);
  }

  //********Parameters********************
  int m_iRotations=atoi(argv[3]);
  double m_fDistThr=atof(argv[4]);

  cout<<"Number of rotations: "<< m_iRotations<<endl;
  cout<<"Distance threshold: "<< m_fDistThr <<endl;

  Molecule<Atom> molModel,molTarget;

  ifstream fileModel(argv[2],ios::in);
  ifstream fileTarget(argv[1],ios::in);

  if(!fileModel){
    cout<<"File "<<argv[1]<< "does not exist."<<endl;
    exit(0);
  }
  if(!fileTarget){
    cout<<"File "<<argv[2]<< "does not exist."<<endl;
    exit(0);
  }

  molModel.readPDBfile(fileModel, PDB::CAlphaSelector());
  molTarget.readPDBfile(fileTarget, PDB::CAlphaSelector());

  //calculate center of mass
  Vector3 vectModelMass(0,0,0);
  for(unsigned int i=0;i<molModel.size();i++) {
    vectModelMass+=molModel[i].position();
  }
  vectModelMass/=molModel.size();

  Vector3 vectTargetMass(0,0,0);
  for(unsigned int i=0;i<molTarget.size();i++) {
    vectTargetMass+=molTarget[i].position();
  }
  vectTargetMass/=molTarget.size();

  //transform the molecules to the center of the coordinate system
  molModel+=(-vectModelMass);
  molTarget+=(-vectTargetMass);

  //insert into hash
  GeomHash <Vector3,int> gHash(3,m_fDistThr);

  //loop over all atoms of target molecule
  for(unsigned int i=0;i<molTarget.size();i++) {
    gHash.insert(molTarget[i].position(),i);
  }

  //choose the best alignment from random rotations
  unsigned int iMaxSize=0;
  RigidTrans3 rtransBest;

  for(int iRand=0; iRand < m_iRotations; iRand++) {
    Matrix3 rotation((drand48()-0.5)*2*3.1415,
                     (drand48()-0.5)*2*3.1415,
                     (drand48()-0.5)*2*3.1415); //choose random rotation

    Match match;
    //adding the pairs of atom that are close enough
    for(unsigned int i=0;i< molModel.size();i++) {
      Vector3 mol_atom = rotation*molModel[i].position();

      HashResult<int> result;
      gHash.query(mol_atom, m_fDistThr, result);

      //check if the atoms inside distance threshold
      HashResult<int>::iterator x_end=result.end(), x;
      for(x = result.begin(); x != x_end; x++) {
        if(mol_atom.dist(molTarget[*x].position()) <= m_fDistThr) {
          float score = (1 / (1+(molModel[i].position() | molTarget[ *x ])));
          match.add( *x , i, score, score );
        }
      }
      result.clear();
    }

    //calculates transformation that is a little better than "rotation"
    match.calculateBestFit(molTarget, molModel);

    if(iMaxSize < match.size() ){
      iMaxSize = match.size();
      rtransBest=match.rigidTrans();
    }
  }

  cout<< "Max Alignment Size: " << iMaxSize <<endl;
  cout<< "Rigid Trans: "<<
    RigidTrans3(Vector3(0,0,0),vectTargetMass)*
    rtransBest*
    RigidTrans3(Vector3(0,0,0),(-vectModelMass))<<endl;

  cout<<overallTime<<endl;
}
