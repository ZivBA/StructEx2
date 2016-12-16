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
    double epsilon = atof(argv[1]);

    cout << "Distance threshold: " << epsilon << endl;

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


    Molecule<Atom> modelAlpha, targetAlpha;
    molModel.select(Molecule<Atom>::CaSelector(), modelAlpha);
    molTarget.select(Molecule<Atom>::CaSelector(), targetAlpha);


    unsigned long scoreY = modelAlpha.size() + 1;
    unsigned long scoreX = targetAlpha.size() + 1;
    // scoring matrix for sequence alignment
    double scoreMatrix[scoreY][scoreX];
    // backstep matrix
    int backMatrix[scoreY][scoreX];

    for (int l = 0; l < scoreY; l++)
    {
        scoreMatrix[l][0] = 0.0;
        backMatrix[l][0] = 0;
    }

    for (int l = 0; l < scoreX; l++)
    {
        scoreMatrix[0][l] = 0.0;
        backMatrix[0][l] = 0;
    }


    RigidTrans3 bestTrans;

    int bestAlignSize = 0;

    for (unsigned int i = 0; i < backBoneModel.size(); i += 3)
    {
        Triangle trigModel = Triangle(backBoneModel[i].position(), backBoneModel[i + 1].position(),
                                      backBoneModel[i + 2].position());

        if (true){
            cout << "did "<< i/3 << " trigModel iterations " << endl;
        }


        for (unsigned int j = 0; j < backBoneTarget.size(); j += 3)
        {

            if (j % 100 ==0){
                cout << "did "<< j/3 << " trigTarget iterations " << endl;
            }

            Triangle trigTarget = Triangle(backBoneTarget[j].position(),
                                           backBoneTarget[j + 1].position(),
                                           backBoneTarget[j + 2].position());

            RigidTrans3 curTrans = trigModel | trigTarget;
            molTarget.select(Molecule<Atom>::BackboneSelector(), tempTarget);


            molModel.select(Molecule<Atom>::CaSelector(), modelAlpha);
            molTarget.select(Molecule<Atom>::CaSelector(), targetAlpha);

            targetAlpha *= curTrans;





            for (unsigned int k = 1; k < scoreY; k++)
            {

                for (unsigned int l = 1; l < scoreX; l++)
                {


                    float dist = modelAlpha[k - 1] | targetAlpha[l - 1];
                    double match = dist < epsilon ? (epsilon - dist) : 0;

                    int traceScore = scoreMatrix[k][l] == match + scoreMatrix[k - 1][l - 1] ? 1 : 0;

                    double score = max3(scoreMatrix[k - 1][l], scoreMatrix[k][l - 1],
                                       match + scoreMatrix[k - 1][l - 1]);

                    if (score == scoreMatrix[k - 1][l]){
                        scoreMatrix[k][l] = score;
                        backMatrix[k][l] = backMatrix[k-1][l];
                    } else if (score == scoreMatrix[k][l-1]) {
                        scoreMatrix[k][l] = score;
                        backMatrix[k][l] = backMatrix[k][l-1];
                    } else{
                        scoreMatrix[k][l] = score;
                        backMatrix[k][l] = backMatrix[k-1][l-1]+traceScore;
                    }

                }
            }

            if (backMatrix[scoreY-1][scoreX-1] > bestAlignSize){
                bestTrans = curTrans;
                bestAlignSize = backMatrix[scoreY-1][scoreX-1];

            }


        }

    }


    molTarget *= bestTrans;

    ofstream myfile;
    myfile.open("transformed.pdb");
    myfile << molTarget;
    myfile.close();

    cout << overallTime << endl;
}
