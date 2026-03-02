/**
 * @file		aux.test.c
 * @brief		Unit tests for aux.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 */

#include "core.h"
#include "test.h"
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>

static unsigned int nextRandU32(unsigned int *state){
	unsigned int x = *state;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	*state = x;
	return x;
}

static unsigned int randomSeed(){
	unsigned int seed = (unsigned int)time(NULL) ^ (unsigned int)clock();
	if(seed == 0){
		seed = 0xA341316Cu;
	}
	return seed;
}

static int randIntInRange(unsigned int *state, int min, int max){
	unsigned int span = (unsigned int)(max - min + 1);
	return min + (int)(nextRandU32(state) % span);
}

static double randDoubleInRange(unsigned int *state, double min, double max){
	double unit = (double)nextRandU32(state) / (double)UINT_MAX;
	return min + (max - min)*unit;
}

static int nearlyEqual(double a, double b, double absTol, double relTol){
	double diff = fabs(a-b);
	double scale = fmax(fabs(a),fabs(b));
	return diff <= absTol + relTol*scale;
}

static int assertIntArrayEq(const int *a, const int *b, long int n, const char *label){
	for(long int i=0;i<n;i++){
		utAssert(a[i]==b[i],"%s mismatch at i=%li: got=%i expected=%i",label,i,a[i],b[i]);
	}
	return 0;
}

static int assertLongArrayEq(const long int *a, const long int *b, long int n, const char *label){
	for(long int i=0;i<n;i++){
		utAssert(a[i]==b[i],"%s mismatch at i=%li: got=%li expected=%li",label,i,a[i],b[i]);
	}
	return 0;
}

static int assertDoubleArrayEq(const double *a, const double *b, long int n, double tol, const char *label){
	for(long int i=0;i<n;i++){
		utAssert(fabs(a[i]-b[i])<tol,"%s mismatch at i=%li: got=%f expected=%f",label,i,a[i],b[i]);
	}
	return 0;
}

static int testStrCatAlloc(){
	char *s = strCatAlloc(4,"PIN","C","-","test");
	utAssert(strcmp(s,"PINC-test")==0,"strCatAlloc returned '%s'",s);
	free(s);
	return 0;
}

static int testAiProd(){
	int arr[] = {2,34,9,1,6,6};
	long int prod = aiProd(arr,6);
	utAssert(prod==22032,"aiProd doesn't work");
	return 0;
}

static int testAEq(){
	int a[] = {2,3,4,5,6};
	int b[] = {2,3,4,5,6};
	long int al[] = {2,3,4,5,6};
	long int bl[] = {2,3,4,5,6};

	utAssert(aiEq(a,b,5),"aiEq is broken");
	utAssert(alEq(al,bl,5),"alEq is broken");

	b[4] = 3;
	bl[4] = 3;

	utAssert(!aiEq(a,b,5),"aiEq is broken");
	utAssert(!alEq(al,bl,5),"alEq is broken");

	return 0;
}

static int testBasicIntVectorOps(){
	const long int n = 5;
	int a[] = {1,2,3,4,5};
	int b[] = {5,4,3,2,1};
	int res[] = {0,0,0,0,0};

	aiAdd(a,b,res,n);
	int addExp[] = {6,6,6,6,6};
	assertIntArrayEq(res,addExp,n,"aiAdd");

	aiMul(a,b,res,n);
	int mulExp[] = {5,8,9,8,5};
	assertIntArrayEq(res,mulExp,n,"aiMul");

	aiScale(a,n,3);
	int scaleExp[] = {3,6,9,12,15};
	assertIntArrayEq(a,scaleExp,n,"aiScale");

	aiShift(a,n,-2);
	int shiftExp[] = {1,4,7,10,13};
	assertIntArrayEq(a,shiftExp,n,"aiShift");

	aiSetAll(a,n,9);
	int setAllExp[] = {9,9,9,9,9};
	assertIntArrayEq(a,setAllExp,n,"aiSetAll");

	aiSet(a,n,7,6,5,4,3);
	int setExp[] = {7,6,5,4,3};
	assertIntArrayEq(a,setExp,n,"aiSet");

	return 0;
}

static int testReductionsAndExtrema(){
	int a[] = {-5,2,7,-3,1};
	long int al[] = {-5,2,7,-3,1};
	double ad[] = {-5.0,2.0,7.0,-3.0,1.0};
	const long int n = 5;

	utAssert(aiMin(a,n)==-5,"aiMin failed");
	utAssert(aiMax(a,n)==7,"aiMax failed");
	utAssert(aiExt(a,n)==7,"aiExt failed");
	utAssert(aiSum(a,n)==2,"aiSum failed");
	utAssert(fabs(aiAvg(a,n)-0.4)<1e-15,"aiAvg failed");

	utAssert(alMin(al,n)==-5,"alMin failed");
	utAssert(alMax(al,n)==7,"alMax failed");
	utAssert(alExt(al,n)==7,"alExt failed");
	utAssert(alSum(al,n)==2,"alSum failed");
	utAssert(fabs(alAvg(al,n)-0.4)<1e-15,"alAvg failed");

	utAssert(fabs(adMin(ad,n)+5.0)<1e-15,"adMin failed");
	utAssert(fabs(adMax(ad,n)-7.0)<1e-15,"adMax failed");
	utAssert(fabs(adExt(ad,n)-7.0)<1e-15,"adExt failed");
	utAssert(fabs(adSum(ad,n)-2.0)<1e-15,"adSum failed");
	utAssert(fabs(adAvg(ad,n)-0.4)<1e-15,"adAvg failed");
	utAssert(fabs(adProd(ad,n)-210.0)<1e-15,"adProd failed");

	return 0;
}

static int testCumOps(){
	int a[] = {5,4,3};
	const long int n = 3;
	int aiProdOut[4];
	int aiSumOut[4];
	long int ailProdOut[4];
	long int ailSumOut[4];

	aiCumProd(a,aiProdOut,n);
	int aiProdExp[] = {1,5,20,60};
	assertIntArrayEq(aiProdOut,aiProdExp,4,"aiCumProd");

	aiCumSum(a,aiSumOut,n);
	int aiSumExp[] = {0,5,9,12};
	assertIntArrayEq(aiSumOut,aiSumExp,4,"aiCumSum");

	ailCumProd(a,ailProdOut,n);
	long int ailProdExp[] = {1,5,20,60};
	assertLongArrayEq(ailProdOut,ailProdExp,4,"ailCumProd");

	ailCumSum(a,ailSumOut,n);
	long int ailSumExp[] = {0,5,9,12};
	assertLongArrayEq(ailSumOut,ailSumExp,4,"ailCumSum");

	return 0;
}

static int testDotProdAndDoubleEq(){
	double a[] = {1.0,2.0,3.0};
	double b[] = {4.0,5.0,6.0};
	double c[] = {1.0,2.0,3.0000000001};
	const long int n = 3;

	utAssert(adDotProd(a,b,n)==32,"adDotProd failed");
	utAssert(adEq(a,c,n,1e-9),"adEq should pass with loose tolerance");
	utAssert(!adEq(a,c,n,1e-12),"adEq should fail with strict tolerance");

	return 0;
}

static int testLongVariantOps(){
	const long int n = 4;
	long int a[] = {2,3,4,5};
	long int b[] = {6,7,8,9};
	long int res[] = {0,0,0,0};

	alAdd(a,b,res,n);
	long int addExp[] = {8,10,12,14};
	assertLongArrayEq(res,addExp,n,"alAdd");

	alMul(a,b,res,n);
	long int mulExp[] = {12,21,32,45};
	assertLongArrayEq(res,mulExp,n,"alMul");

	alScale(a,n,2);
	long int scaleExp[] = {4,6,8,10};
	assertLongArrayEq(a,scaleExp,n,"alScale");

	alSet(a,n,9L,8L,7L,6L);
	long int setExp[] = {9,8,7,6};
	assertLongArrayEq(a,setExp,n,"alSet");

	return 0;
}

static int testDoubleVectorOps(){
	const long int n = 3;
	double a[] = {1.0,-2.0,0.5};
	double b[] = {2.0,3.0,4.0};
	double res[] = {0.0,0.0,0.0};

	adAdd(a,b,res,n);
	double addExp[] = {3.0,1.0,4.5};
	assertDoubleArrayEq(res,addExp,n,1e-15,"adAdd");

	adMul(a,b,res,n);
	double mulExp[] = {2.0,-6.0,2.0};
	assertDoubleArrayEq(res,mulExp,n,1e-15,"adMul");

	adScale(a,n,2.0);
	double scaleExp[] = {2.0,-4.0,1.0};
	assertDoubleArrayEq(a,scaleExp,n,1e-15,"adScale");

	adShift(a,n,1.5);
	double shiftExp[] = {3.5,-2.5,2.5};
	assertDoubleArrayEq(a,shiftExp,n,1e-15,"adShift");

	adSet(a,n,7.0,8.0,9.0);
	double setExp[] = {7.0,8.0,9.0};
	assertDoubleArrayEq(a,setExp,n,1e-15,"adSet");

	return 0;
}

static int testRandomAiAddProperties(){
	const int nCases = 300;
	unsigned int seed = randomSeed();
	unsigned int rng = seed;

	for(int c=0;c<nCases;c++){
		long int n = randIntInRange(&rng,1,64);

		int *a = malloc((size_t)n*sizeof(*a));
		int *b = malloc((size_t)n*sizeof(*b));
		int *ab = malloc((size_t)n*sizeof(*ab));
		int *ba = malloc((size_t)n*sizeof(*ba));
		utAssert(a && b && ab && ba,"malloc failed in random aiAdd test; seed=%u case=%i",seed,c);

		for(long int i=0;i<n;i++){
			a[i] = randIntInRange(&rng,-1000,1000);
			b[i] = randIntInRange(&rng,-1000,1000);
		}

		aiAdd(a,b,ab,n);
		aiAdd(b,a,ba,n);

		for(long int i=0;i<n;i++){
			utAssert(ab[i]==ba[i],
				"aiAdd commutativity failed; seed=%u case=%i n=%li i=%li a=%i b=%i ab=%i ba=%i",
				seed,c,n,i,a[i],b[i],ab[i],ba[i]);
		}

		{
			long int sumAB = aiSum(ab,n);
			long int sumA = aiSum(a,n);
			long int sumB = aiSum(b,n);
			utAssert(sumAB==sumA+sumB,
				"aiAdd/sum consistency failed; seed=%u case=%i n=%li sumAB=%li sumA=%li sumB=%li",
				seed,c,n,sumAB,sumA,sumB);
		}

		free(a);
		free(b);
		free(ab);
		free(ba);
	}

	return 0;
}

static int testRandomAdAddCommutativity(){
	const int nCases = 300;
	unsigned int seed = randomSeed();
	unsigned int rng = seed;
	const double absTol = 1e-12;
	const double relTol = 1e-12;

	for(int c=0;c<nCases;c++){
		long int n = randIntInRange(&rng,1,64);

		double *a = malloc((size_t)n*sizeof(*a));
		double *b = malloc((size_t)n*sizeof(*b));
		double *ab = malloc((size_t)n*sizeof(*ab));
		double *ba = malloc((size_t)n*sizeof(*ba));
		utAssert(a && b && ab && ba,"malloc failed in random adAdd test; seed=%u case=%i",seed,c);

		for(long int i=0;i<n;i++){
			a[i] = randDoubleInRange(&rng,-1e6,1e6);
			b[i] = randDoubleInRange(&rng,-1e6,1e6);
		}

		adAdd(a,b,ab,n);
		adAdd(b,a,ba,n);

		for(long int i=0;i<n;i++){
			utAssert(nearlyEqual(ab[i],ba[i],absTol,relTol),
				"adAdd commutativity failed; seed=%u case=%i n=%li i=%li a=%g b=%g ab=%g ba=%g",
				seed,c,n,i,a[i],b[i],ab[i],ba[i]);
		}

		free(a);
		free(b);
		free(ab);
		free(ba);
	}

	return 0;
}

// All tests for aux.c is contained in this function
void testAux(){
	utRun(&testStrCatAlloc);
	utRun(&testAiProd);
	utRun(&testAEq);
	utRun(&testBasicIntVectorOps);
	utRun(&testReductionsAndExtrema);
	utRun(&testCumOps);
	utRun(&testDotProdAndDoubleEq);
	utRun(&testLongVariantOps);
	utRun(&testDoubleVectorOps);
	utRun(&testRandomAiAddProperties);
	utRun(&testRandomAdAddCommutativity);
}
