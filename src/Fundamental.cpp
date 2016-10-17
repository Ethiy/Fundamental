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
#include <limits>

// STD C includes
#include <cstdlib>
#include <ctime>
#include <cmath>

// Namespaces
using namespace Imagine;
using namespace std;


// Constants
static const double BETA = 0.01; // Probability of failure
static const double NORM = 1000; // Normalization constant


// Matches structure
struct Match {
    double x1, y1, x2, y2;
};

// Functions prototypes
void algoSIFT(Image<Color,2>, Image<Color,2>, vector<Match>&);

int random_index(int);
int _sample( vector<int>&);
vector<int> sampler( int, int);
bool is_square(size_t, size_t&);
Matrix<double> to_matrix(Vector<double>);
double epipolar_distance(Match, Matrix<double>);
vector<int> infere_inliers( int, vector<Match>&, vector<int>, double);
void estimateNiter(int&, int, int, int);
Matrix<double> estimateF(vector<Match>);
Matrix<double> computeF(vector<Match>&); 
void displayEpipolar(Image<Color>, Image<Color>, Matrix<double>);
void line( Vector<double>, IntPoint2&, IntPoint2&, int, int, int);


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

bool is_square(size_t n, size_t& m)
{
    m = static_cast<size_t>( sqrt( static_cast<double_t>(n)) );
    return abs(n - pow( m, 2)) < numeric_limits<double>::epsilon();
}

Matrix<double> to_matrix(Vector<double> f)
{
    size_t n = 0;
    assert(is_square( f.size(), n)); // making sure f can be represented as a Matrix

    Matrix<double> F(n,n);
    for(size_t i=0; i<n; i++)
        for(size_t j=0; j<n; j++)
            F(i,j) = f[3*i + j];
    return F;
}

double epipolar_distance(Match m, Matrix<double> F)
{
    Vector<double> u(3); u[0] = m.x1; u[1] = m.y1; u[2] = 1;
    Vector<double> v(3); v[0] = m.x2; v[1] = m.y2; v[2] = 1;
    v = transpose(F) * v;
    v = v / ( sqrt( pow(v[0],2.) + pow( v[1],2.)));
    return abs(u*v);
}

vector<int> infere_inliers( int n_samples, vector<Match>& matches, vector<int> sampled_matches, double sigma)
{
    assert( n_samples == static_cast<int>( sampled_matches.size()));
    vector<int> inliers;

    // Defining System
    Matrix<double> A = Matrix<double>::Zero(n_samples + 1,9);
    Matrix<double> F;

    for(int i = 0; i < n_samples; i++)
    {
        A(i,0) = matches[ sampled_matches[i]].x1 * matches[ sampled_matches[i]].x2;
        A(i,1) = matches[ sampled_matches[i]].y1 * matches[ sampled_matches[i]].x2;
        A(i,2) = matches[ sampled_matches[i]].x2;
        A(i,3) = matches[ sampled_matches[i]].x1 * matches[ sampled_matches[i]].y2;
        A(i,4) = matches[ sampled_matches[i]].y1 * matches[ sampled_matches[i]].y2;
        A(i,5) = matches[ sampled_matches[i]].y2;
        A(i,6) = matches[ sampled_matches[i]].x1;
        A(i,7) = matches[ sampled_matches[i]].y1;
        A(i,8) = 1;
    }

    // To optimize memory we put all auxilary variables in a scope
    {
        Matrix<double> U,V;
        Vector<double> S;
        svd(A,U,S,V);
        V = transpose(V);
        Vector<double> f = static_cast<Vector<double>>( V.getSubMat(0,9,8,1));
        F = to_matrix(f);

        // Enforce rank 2
        svd(F,U,S,V);
        cout << S << endl << flush;
        S[2] = 0;
        F = U * Diagonal(S) * V;
    }

    // Find Inliers
    for( int index = 0; index < static_cast<int>(matches.size()); ++index)
        if( epipolar_distance( matches[index], F) <= sigma)
            inliers.push_back(index);

    return inliers;
}

void estimateNiter(int& Niter, int n_inliers, int n_samples, int n_matches)
{
    assert(n_inliers != 0); // otherwise Niter == inf
    double aux = pow( static_cast<double>(n_inliers)/static_cast<double>(n_matches), static_cast<double>(n_samples));
    if( abs(aux) > numeric_limits<double>::epsilon()) // otherwise Niter == inf
        Niter = max( static_cast<int>( ceil( log(BETA)/log(1 - aux ))), 0);
}

Matrix<double> estimateF( vector<Match> matches)
{
    Matrix<double> F;
    int n = static_cast<int>(matches.size());

    // Defining System
    Matrix<double> A = Matrix<double>::Zero(n,9);
    Vector<double> B(n);
    B.fill(0);

    for(int i = 0; i < n; i++)
    {
        A(i,0) = matches[i].x1 * matches[i].x2;
        A(i,1) = matches[i].y1 * matches[i].x2;
        A(i,2) = matches[i].x2;
        A(i,3) = matches[i].x1 * matches[i].y2;
        A(i,4) = matches[i].y1 * matches[i].y2;
        A(i,5) = matches[i].y2;
        A(i,6) = matches[i].x1;
        A(i,7) = matches[i].y1;
        A(i,8) = 1;
    }
    // Enforce rank 2
    {
        Matrix<double> U,V;
        Vector<double> S;
        svd(A,U,S,V);
        V = transpose(V);
        Vector<double> f = static_cast<Vector<double>>( V.getSubMat(0,9,8,1));
        cout << "   Mean square Error = " << norm(A*f) << endl << flush;
        F = to_matrix(f);
        svd(F,U,S,V);
        S[2] = 0;
        F = U * Diagonal(S) * V;
        S[0] = 1.0/NORM; S[1] = 1.0/NORM; S[2] = 1.0;
        F = Diagonal(S)*F*Diagonal(S);
    }
    return F;
}

void line( Vector<double> v, IntPoint2& left, IntPoint2& right, int w, int h, int view)
{
    assert( view == 0 || view == 1);
    if( abs(v[1]) > numeric_limits<double>::epsilon() )
    {
        left[0] = view * w; right[0] = w * (1 + view);

        if( abs(v[0]) > numeric_limits<double>::epsilon() )
        {
            left[1] = - static_cast<int>( round(v[2]/v[1]));
            if(left[1]<0)
            {
                left[0] = - static_cast<int>( round(v[2] / v[0])) + view * w;
                left[1] = 0;
            }
            right[1] = - static_cast<int>( round(v[0] / v[1])) * w - static_cast<int>( round(v[2] / v[1]));
            if(right[1]>h)
            {
                right[1] = view * w;
                right[0] = - static_cast<int>( round(v[2] / v[0])) - static_cast<int>( round(v[1] / v[0])) * h;
            }
        }
        else
        {
            left[1] = - static_cast<int>( round(v[2] /  v[1]));
            right[1] = - static_cast<int>( round(v[2] / v[1]));
        }
    }
    else
    {
        left[0] = - static_cast<int>( round(v[2] /  v[0])) + view * w;
        left[1] = 0;
        right[0] = - static_cast<int>( round(v[2] /  v[0])) + view * w;
        right[1] = h;
    }
}

// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
Matrix<double> computeF(vector<Match>& matches) 
{
    const double distMax(1.5); // Pixel error for inlier/outlier discrimination
    int Niter(100000); // Adjusted dynamically

    Matrix<double> F;
    vector<int> inliers;

    const int n_samples(8);
    const int n_matches = static_cast<int>( matches.size());

    cout << "Normalizing..." << endl << flush;
    // 0. Normalization
    vector<Match> normalized_matches(matches);
    for(size_t i=0; i<normalized_matches.size(); ++i)
    {
        normalized_matches[i].x1 /= NORM;
        normalized_matches[i].x2 /= NORM;
        normalized_matches[i].y1 /= NORM;
        normalized_matches[i].y2 /= NORM;
    }

    cout << "Ransac..." << endl << flush;
    // Iterating
    int n_inliers = n_samples; // In the worse case scenario there is at least the minimum points to estimate F
    int counter(0);
    vector<int> sampled_matches( n_samples, 0);
    while( counter < Niter )
    {
        cout << "Iteration: " << counter << endl << flush;
        // 1. Sample matches
        cout << "   Sample points..." << endl << flush;
        sampled_matches = sampler( n_samples, n_matches);

        // 2. Infere inliers
        cout << "   Inliers within model..." << endl << flush;
        inliers = infere_inliers( n_samples, normalized_matches, sampled_matches, distMax);
        // 3. Update Niter
        if( static_cast<int>( inliers.size()) >= n_inliers )
            n_inliers = static_cast<int>( inliers.size());
        cout << "   Number of inliers : " << n_inliers << endl << flush;

        estimateNiter( Niter, n_inliers, n_samples, n_matches);
        cout << "   Niter estimate: " << Niter << endl << flush;
        counter++;
    }

    // Keeping inliers only
    cout << "Keeping only inliers..." << endl << flush;
    vector<Match> all(matches);
    matches.clear();
    for(size_t i=0; i<inliers.size(); i++)
        matches.push_back(all[inliers[i]]);

    normalized_matches.clear();
    for(size_t i=0; i<inliers.size(); i++)
        normalized_matches.push_back(all[inliers[i]]);

    // Estimating F
    cout << "Estimating F ..." << endl << flush;
    F = estimateF(normalized_matches);
    return F;
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2, Matrix<double> F) 
{
    while(true) 
    {
        Vector<double> v(3);
        v[2] = 1;
        int w = I1.width();
        int x,y;
        if(getMouse(x,y) == 3)
            break;
        cout << "Point clicked = " << x << ' ' << y << endl << flush;
        int view = static_cast<int>(x<w);
        cout << "The epipolar view is in Image : " << view + 1 << endl << flush;
        v[1] = y;
        v[0] = x;
        v = ( static_cast<double>(view) * transpose(F) + static_cast<double>(1 - view) * F) * v;
        cout << v << endl;
        IntPoint2 left,right;
        int h = (view * I2.height() + (1-view) * I1.height());
        line( v, left, right, w, h, view);
        cout << left << " and  " << right << endl;
        drawLine( left, right, RED, 2); 
    }
}

int main(int argc, char** argv)
{
    cout << "Reading Images..." << endl << flush;
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

    cout << "Images displayed." << endl << flush;

    vector<Match> matches;
    cout << "Running SIFT Algorithm..." << endl << flush;
    algoSIFT(I1, I2, matches);
    cout << " matches: " << matches.size() << endl;
    cout << "Waiting for a click..." << endl << flush;
    click();
    
    cout << "Computing F..." << endl << flush;
    Matrix<double> F = computeF(matches);
    cout << "F="<< endl << F << flush;

    // Redisplay with matches
    cout << "Redisplaying with matches..." << endl << flush;
    display(I1,0,0);
    display(I2,w,0);
    for(size_t i=0; i<matches.size(); i++) {
        Color c( static_cast<unsigned char>( random_index(256)), static_cast<unsigned char>(random_index(256)), static_cast<unsigned char>(random_index(256)));
        fillCircle( static_cast<int>(matches[i].x1), static_cast<int>(matches[i].y1), 2, c);
        fillCircle(static_cast<int>(matches[i].x2+w), static_cast<int>(matches[i].y2), 2, c);        
    }
    cout << "Waiting for a click..." << endl << flush;
    click();

    // Redisplay without SIFT points
    display(I1,0,0);
    display(I2,w,0);
    cout << "Displaying epipolar lines..." << endl << flush;
    displayEpipolar(I1, I2, F);

    endGraphics();
    return EXIT_SUCCESS;
}
