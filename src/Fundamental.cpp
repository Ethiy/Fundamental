// Imagine++ project
// Project:  Fundamental
// Author:   Pascal Monasse
// Modified by: Oussama Ennafii
// Date:     2016/10/12

// local includes
#include "../Imagine/Features.h"

// Imagine includes
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>

// STD C++ includes
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>

// STD C includes
#include <cstdlib>
#include <ctime>
#include <cmath>

// Namespaces
using namespace Imagine;
using namespace std;


// Constants
static const float BETA = 0.01f; // Probability of failure

// Structures
struct Match {
    float x1, y1, x2, y2;
};

// Functions

// Display SIFT points and fill vector of point correspondences
void algoSIFT(Image<Color,2> I1, Image<Color,2> I2, vector<Match>& matches) 
{
    // Find interest points
    SIFTDetector D;
    D.setFirstOctave(-1);
    Array<SIFTDetector::Feature> feats1 = D.run(I1);
    drawFeatures(feats1, Coords<2>(0,0));
    cout << "Im1: " << feats1.size() << flush;
    Array<SIFTDetector::Feature> feats2 = D.run(I2);
    drawFeatures(feats2, Coords<2>(I1.width(),0));
    cout << " Im2: " << feats2.size() << flush;

    const double MAX_DISTANCE = 100.0*100.0;
    for(size_t i=0; i < feats1.size(); i++) {
        SIFTDetector::Feature f1=feats1[i];
        for(size_t j=0; j < feats2.size(); j++) {
            double d = squaredDist(f1.desc, feats2[j].desc);
            if(d < MAX_DISTANCE) {
                Match m;
                m.x1 = f1.pos.x();
                m.y1 = f1.pos.y();
                m.x2 = feats2[j].pos.x();
                m.y2 = feats2[j].pos.y();
                matches.push_back(m);
            }
        }
    }
}

int random_index(int len)
{
    random_device device;
    mt19937 gen(device());
    uniform_int_distribution<int> X( 0, len-1);
    return X(gen);
}

int _sample( vector<int>& v)
{
    assert(v.size());
    int k = random_index( static_cast<int>( v.size()) );
    int value = v[k];
    v.erase( v.begin() + k);
    return value;
}

vector<int> sampler( int n_samples, int n_matches)
{
    assert( n_samples <= n_matches && n_samples );

    vector<int> samples( n_samples, 0);
    vector<int> v(n_matches);
    iota(v.begin(), v.end(), 0);
    
    for(int iterator = 0; iterator < n_samples; iterator++)
        samples[iterator] = _sample(v);
    return samples;
}

vector<int> infereF( vector<int> sampled_matches, float sigma, FMatrix<float,3,3> bestF)
{
    vector<int> inliers;
    return inliers;
}

void estimateNiter(int& Niter, int n_inliers, int n_samples, int n_matches)
{
    assert(n_inliers != 0); // otherwise Niter == inf
    Niter = max( static_cast<int>( ceil( log(BETA)/log(1- pow( static_cast<float>(n_inliers)/static_cast<float>(n_matches), static_cast<float>(n_samples))))), 0);
}

// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
FMatrix<float,3,3> computeF(vector<Match>& matches) 
{
    const float distMax(1.5f); // Pixel error for inlier/outlier discrimination
    int Niter(100000); // Adjusted dynamically
    FMatrix<float,3,3> bestF;
    vector<int> bestInliers;
    // --------------- TODO ------------
    // DO NOT FORGET NORMALIZATION OF points
    const int n_samples(8);
    const int n_matches = static_cast<int>( matches.size());
    int n_inliers = n_samples; // In the worse case scenario there is at least the minimum points to estimate F

    // Iterating
    int counter(0);
    vector<int> sampled_matches( n_samples, 0);
    while( counter < Niter )
    {
        sampled_matches = sampler( n_samples, n_matches);
        bestInliers = infereF( sampled_matches, distMax, bestF);
        if( static_cast<int>( bestInliers.size()) >= n_inliers )
            n_inliers = static_cast<int>( bestInliers.size());
        counter++;
        estimateNiter( Niter, n_inliers, n_samples, n_matches);
    }  


    // Updating matches with inliers only
    vector<Match> all=matches;
    matches.clear();
    for(size_t i=0; i<bestInliers.size(); i++)
        matches.push_back(all[bestInliers[i]]);
    return bestF;
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2,
                     const FMatrix<float,3,3>& F) {
    while(true) {
        int x,y;
        if(getMouse(x,y) == 3)
            break;
        // --------------- TODO ------------
    }
}

int main(int argc, char* argv[])
{
    srand((unsigned int)time(0));

    const char* s1 = argc>1? argv[1]: srcPath("ressources/images/im1.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("ressources/images/im2.jpg");

    // Load and display images
    Image<Color,2> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load images" << endl;
        return EXIT_FAILURE;
    }

    int w = I1.width();
    openWindow(2*w, I1.height());
    display(I1,0,0);
    display(I2,w,0);

    vector<Match> matches;
    algoSIFT(I1, I2, matches);
    cout << " matches: " << matches.size() << endl;
    click();
    
    FMatrix<float,3,3> F = computeF(matches);
    cout << "F="<< endl << F << flush;

    // Redisplay with matches
    display(I1,0,0);
    display(I2,w,0);
    for(size_t i=0; i<matches.size(); i++) {
        Color c( rand()%256, rand()%256, rand()%256);
        fillCircle(matches[i].x1+0, matches[i].y1, 2, c);
        fillCircle(matches[i].x2+w, matches[i].y2, 2, c);        
    }
    click();

    // Redisplay without SIFT points
    display(I1,0,0);
    display(I2,w,0);
    displayEpipolar(I1, I2, F);

    endGraphics();
    return EXIT_SUCCESS;
}
