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
#include <Triangle.h>


int main(int argc, char *argv[])
{
    Timer overallTime("Overall Run Time");
    overallTime.start();

    if (argc != 4)
    {
        cerr << "Usage: " << argv[0] << " epsilon pdb1 pdb2"
             << endl;
        exit(1);
    }

    //********Parameters********************
    const float epsilon = (const float) atof(argv[1]);

//    cout << "Distance threshold: " << epsilon << endl;

    Molecule<Atom> molModel, molTarget, backBoneModel, backBoneTarget, tempTarget;

    ifstream fileModel(argv[2], ios::in);
    ifstream fileTarget(argv[3], ios::in);

    if (!fileModel)
    {
        cout << "File " << argv[2] << "does not exist." << endl;
        exit(0);
    }
    if (!fileTarget)
    {
        cout << "File " << argv[3] << "does not exist." << endl;
        exit(0);
    }

    molModel.readPDBfile(fileModel);
    molTarget.readPDBfile(fileTarget);

    molModel.select(Molecule<Atom>::BackboneSelector(), backBoneModel);
    molTarget.select(Molecule<Atom>::BackboneSelector(), backBoneTarget);


    Molecule<Atom> alphaModel, alphaTarget;
    molModel.select(Molecule<Atom>::CaSelector(), alphaModel);
    molTarget.select(Molecule<Atom>::CaSelector(), alphaTarget);

    GeomHash<Vector3, int> gHash(3, epsilon);

    for (unsigned int i = 0; i < alphaModel.size(); i++)
    {
        gHash.insert(alphaModel[i].position(), i);
    }


//    unsigned long scoreY = alphaModel.size() + 1;
//    unsigned long scoreX = alphaTarget.size() + 1;
//    // scoring matrix for sequence alignment
//    double scoreMatrix[scoreY][scoreX];
//    // backstep matrix
//    int backMatrix[scoreY][scoreX];
//
//    for (int l = 0; l < scoreY; l++)
//    {
//        scoreMatrix[l][0] = 0.0;
//        backMatrix[l][0] = 0;
//    }
//
//    for (int l = 0; l < scoreX; l++)
//    {
//        scoreMatrix[0][l] = 0.0;
//        backMatrix[0][l] = 0;
//    }


    RigidTrans3 bestTrans;

    int bestAlignSize = 0;

    float RMSD = 0.0;



    for (unsigned int i = 0; i < backBoneModel.size(); i += 4)
    {

        Triangle trigModel = Triangle(backBoneModel[i].position(), backBoneModel[i + 1].position(),
                                      backBoneModel[i + 2].position());

//        if (i % 30 == 0)
//        {
//            cout << "did " << i / 3 << " trigModel iterations " << endl;
//        }


        for (unsigned int j = 0; j < backBoneTarget.size(); j += 4)
        {

//            if (j % 360 ==0){
//                cout << "did "<< j/3 << " trigTarget iterations " << endl;
//            }

            Molecule<Atom> target4RMSD;
            Molecule<Atom> model4RMSD;

            Triangle trigTarget = Triangle(backBoneTarget[j].position(),
                                           backBoneTarget[j + 1].position(),
                                           backBoneTarget[j + 2].position());

            RigidTrans3 curTrans = trigModel | trigTarget;

            Match match;
            HashResult<int> result;
            int lastIndex = 0;
            int lastMatchSize = 0;

            for (unsigned int f = 0; f < alphaTarget.size(); f++)
            {
                Vector3 mol_atom = curTrans * alphaTarget[f].position();

                gHash.query(mol_atom, epsilon, result);
                HashResult<int>::iterator x;

                for (x = result.begin(); x != result.end(); x++)
                {
                    if (mol_atom.dist(alphaModel[*x].position()) <= epsilon)
                    {
                        float score = (1 / (1 + (alphaTarget[f].position() | alphaModel[*x])));
                        if (/* *x >= lastIndex */ true)
                        {
                            lastIndex = *x;
                            match.add(*x, f, score, score);
                            model4RMSD.add(alphaModel[*x]);
                            target4RMSD.add(Atom(mol_atom, f));
                        }
                    }
                }
                result.clear();
            }
            //calculates transformation that is a little better than "rotation"

//            match.calculateBestFit(alphaTarget,alphaModel);

            if (bestAlignSize < match.size())
            {
                bestAlignSize = match.size();
                bestTrans = curTrans;
                RMSD = match.calculateRMSD(model4RMSD, target4RMSD);

            }
        }

    }

    Match match;
    molTarget *= bestTrans;
//    match.calculateBestFit(molModel,molTarget);
//    RMSD = match.calculateRMSD(molModel,molTarget);
//    cout << bestTrans << endl;
    ofstream myfile;
    myfile.open("./transformed.pdb");
    myfile << molTarget;
    myfile.close();
    cout << bestAlignSize << " " << RMSD << " " << bestTrans << endl;
    cout << overallTime << endl;
}
