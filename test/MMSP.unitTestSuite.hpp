/* MMSP.unitTestSuite.hpp
** Use CxxTest to build unit tests for MMSP classes
** http://cxxtest.com
*/
#include"MMSP.hpp"
#include<cxxtest/TestSuite.h>

class scalarTestSuite : public CxxTest::TestSuite
{
private:
	MMSP::scalar<char>   sc;
	MMSP::scalar<short>  ss;
	MMSP::scalar<int>    si;
	MMSP::scalar<long>   sl;
	MMSP::scalar<float>  sf;
	MMSP::scalar<double> sd;

public:
	void setUp() {
		sc = 2;
		ss = 2;
		si = 2;
		sl = 2;
		sf = 2.0;
		sd = 2.0;
	}
	void testSize(void) {
		TS_ASSERT_EQUALS(sc.buffer_size(),sizeof(char));
		TS_ASSERT_EQUALS(ss.buffer_size(),sizeof(short));
		TS_ASSERT_EQUALS(si.buffer_size(),sizeof(int));
		TS_ASSERT_EQUALS(sl.buffer_size(),sizeof(long));
		TS_ASSERT_EQUALS(sf.buffer_size(),sizeof(float));
		TS_ASSERT_EQUALS(sd.buffer_size(),sizeof(double));
	}
	void testAddition(void) {
		TS_ASSERT_EQUALS(sc+sc,char(4));
		TS_ASSERT_EQUALS(ss+ss,short(4));
		TS_ASSERT_EQUALS(si+si,int(4));
		TS_ASSERT_EQUALS(sl+sl,long(4));
		TS_ASSERT_EQUALS(sf+sf,float(4));
		TS_ASSERT_EQUALS(sd+sd,double(4));
	}
	void testSquares(void) {
		TS_ASSERT_EQUALS(sc*sc,char(4));
		TS_ASSERT_EQUALS(ss*ss,short(4));
		TS_ASSERT_EQUALS(si*si,int(4));
		TS_ASSERT_EQUALS(sl*sl,long(4));
		TS_ASSERT_EQUALS(sf*sf,float(4));
		TS_ASSERT_EQUALS(sd*sd,double(4));
	}
	void testMultiplication(void) {
		TS_ASSERT_EQUALS(2*sc,char(4));
		TS_ASSERT_EQUALS(2*ss,short(4));
		TS_ASSERT_EQUALS(2*si,int(4));
		TS_ASSERT_EQUALS(2*sl,long(4));
		TS_ASSERT_EQUALS(2*sf,float(4));
		TS_ASSERT_EQUALS(2*sd,double(4));
	}
};

class vectorTestSuite : public CxxTest::TestSuite
{
private:
	MMSP::vector<char>   vc;
	MMSP::vector<short>  vs;
	MMSP::vector<int>    vi;
	MMSP::vector<long>   vl;
	MMSP::vector<float>  vf;
	MMSP::vector<double> vd;

public:
	void setUp() {
		MMSP::vector<int> temp(3,2);
		vc = temp;
		vs = temp;
		vi = temp;
		vl = temp;
		vf = temp;
		vd = temp;
	}
	void testValue(void) {
		for (int i=0; i<3; i++) {
			TS_ASSERT_EQUALS(vc[i],char(2));
			TS_ASSERT_EQUALS(vs[i],short(2));
			TS_ASSERT_EQUALS(vi[i],int(2));
			TS_ASSERT_EQUALS(vl[i],long(2));
			TS_ASSERT_EQUALS(vf[i],float(2));
			TS_ASSERT_EQUALS(vd[i],double(2));
		}
	}
	void testSize(void) {
		TS_ASSERT_EQUALS(vc.length(),3);
		TS_ASSERT_EQUALS(vs.length(),3);
		TS_ASSERT_EQUALS(vi.length(),3);
		TS_ASSERT_EQUALS(vl.length(),3);
		TS_ASSERT_EQUALS(vf.length(),3);
		TS_ASSERT_EQUALS(vd.length(),3);

		TS_ASSERT_EQUALS(vc.buffer_size(),sizeof(int)+3*sizeof(char));
		TS_ASSERT_EQUALS(vs.buffer_size(),sizeof(int)+3*sizeof(short));
		TS_ASSERT_EQUALS(vi.buffer_size(),4*sizeof(int));
		TS_ASSERT_EQUALS(vl.buffer_size(),sizeof(int)+3*sizeof(long));
		TS_ASSERT_EQUALS(vf.buffer_size(),sizeof(int)+3*sizeof(float));
		TS_ASSERT_EQUALS(vd.buffer_size(),sizeof(int)+3*sizeof(double));

		TS_ASSERT_THROWS_ANYTHING(vc[4]);
		TS_ASSERT_THROWS_ANYTHING(vs[4]);
		TS_ASSERT_THROWS_ANYTHING(vi[4]);
		TS_ASSERT_THROWS_ANYTHING(vl[4]);
		TS_ASSERT_THROWS_ANYTHING(vf[4]);
		TS_ASSERT_THROWS_ANYTHING(vd[4]);
	}
	void testAddition(void) {
		MMSP::vector<int> temp(3,4);
		for (int i=0; i<3; i++) {
			TS_ASSERT_EQUALS((vc+vc)[i],temp[i]);
			TS_ASSERT_EQUALS((vs+vs)[i],temp[i]);
			TS_ASSERT_EQUALS((vi+vi)[i],temp[i]);
			TS_ASSERT_EQUALS((vl+vl)[i],temp[i]);
			TS_ASSERT_EQUALS((vf+vf)[i],temp[i]);
			TS_ASSERT_EQUALS((vd+vd)[i],temp[i]);
		}
	}
	void testSquares(void) {
		// v*v is the inner (dot) product
		MMSP::vector<int> temp(3,4);
		TS_ASSERT_EQUALS(vc*vc,12);
		TS_ASSERT_EQUALS(vs*vs,12);
		TS_ASSERT_EQUALS(vi*vi,12);
		TS_ASSERT_EQUALS(vl*vl,12);
		TS_ASSERT_EQUALS(vf*vf,12);
		TS_ASSERT_EQUALS(vd*vd,12);
	}
	void testMultiplication(void) {
		MMSP::vector<int> temp(3,4);
		for (int i=0; i<3; i++) {
			TS_ASSERT_EQUALS((2*vc)[i],temp[i]);
			TS_ASSERT_EQUALS((2*vi)[i],temp[i]);
			TS_ASSERT_EQUALS((2*vl)[i],temp[i]);
			TS_ASSERT_EQUALS((2*vf)[i],temp[i]);
			TS_ASSERT_EQUALS((2*vd)[i],temp[i]);
		}
	}
};

class sparseTestSuite : public CxxTest::TestSuite
{
private:
	MMSP::sparse<char>   svc;
	MMSP::sparse<short>  svs;
	MMSP::sparse<int>    svi;
	MMSP::sparse<long>   svl;
	MMSP::sparse<float>  svf;
	MMSP::sparse<double> svd;
	template<typename T> void sparseInit(MMSP::sparse<T>& s) {
		MMSP::sparse<T> temp;
		for (int i=0; i<3; i++)
			temp.set(2*i) = 2;
		s = temp;
	}

public:
	void setUp() {
		sparseInit(svc);
		sparseInit(svs);
		sparseInit(svi);
		sparseInit(svl);
		sparseInit(svf);
		sparseInit(svd);
	}
	void testValue(void) {
		for (int i=0; i<6; i+=2) {
			TS_ASSERT_EQUALS(svc[i],char(2));
			TS_ASSERT_EQUALS(svs[i],short(2));
			TS_ASSERT_EQUALS(svi[i],int(2));
			TS_ASSERT_EQUALS(svl[i],long(2));
			TS_ASSERT_EQUALS(svf[i],float(2));
			TS_ASSERT_EQUALS(svd[i],double(2));
		}
		for (int i=1; i<5; i+=2) {
			TS_ASSERT_EQUALS(svc[i],char(0));
			TS_ASSERT_EQUALS(svs[i],short(0));
			TS_ASSERT_EQUALS(svi[i],int(0));
			TS_ASSERT_EQUALS(svl[i],long(0));
			TS_ASSERT_EQUALS(svf[i],float(0));
			TS_ASSERT_EQUALS(svd[i],double(0));
		}
	}
	void testSize(void) {
		TS_ASSERT_EQUALS(svc.length(),3);
		TS_ASSERT_EQUALS(svs.length(),3);
		TS_ASSERT_EQUALS(svi.length(),3);
		TS_ASSERT_EQUALS(svl.length(),3);
		TS_ASSERT_EQUALS(svf.length(),3);
		TS_ASSERT_EQUALS(svd.length(),3);

		TS_ASSERT_DIFFERS(sizeof(int)+sizeof(char),sizeof(MMSP::item<char>)); // that is unexpected
		TS_ASSERT_EQUALS(sizeof(int)+sizeof(int),sizeof(MMSP::item<char>)); // that is unexpected
		TS_ASSERT_EQUALS(svc.buffer_size(),sizeof(int)+3*sizeof(int)+3*sizeof(char)+3*(sizeof(int)-sizeof(char)));
		TS_ASSERT_EQUALS(svs.buffer_size(),sizeof(int)+3*sizeof(int)+3*sizeof(short)+3*(sizeof(int)-sizeof(short)));
		TS_ASSERT_EQUALS(svi.buffer_size(),sizeof(int)+3*sizeof(int)+3*sizeof(int));
		TS_ASSERT_EQUALS(svl.buffer_size(),sizeof(int)+3*sizeof(int)+3*sizeof(long)+3*(sizeof(long)-sizeof(int)));
		TS_ASSERT_EQUALS(svf.buffer_size(),sizeof(int)+3*sizeof(int)+3*sizeof(float));
		TS_ASSERT_EQUALS(svd.buffer_size(),sizeof(int)+3*sizeof(int)+3*sizeof(double)+3*(sizeof(double)-sizeof(int)));
		// More compactly:
		TS_ASSERT_EQUALS(svc.buffer_size(),sizeof(int)+3*sizeof(int)+3*sizeof(int));
		TS_ASSERT_EQUALS(svs.buffer_size(),sizeof(int)+3*sizeof(int)+3*sizeof(int));
		TS_ASSERT_EQUALS(svi.buffer_size(),sizeof(int)+3*sizeof(int)+3*sizeof(int));
		TS_ASSERT_EQUALS(svl.buffer_size(),sizeof(int)+3*sizeof(long)+3*sizeof(long));
		TS_ASSERT_EQUALS(svf.buffer_size(),sizeof(int)+3*sizeof(int)+3*sizeof(float));
		TS_ASSERT_EQUALS(svd.buffer_size(),sizeof(int)+3*sizeof(double)+3*sizeof(double));
	}
	void testAddition(void) {
		MMSP::sparse<int> temp;
		for (int i=0; i<3; i++)
			temp.set(2*i) = 4;
		TS_ASSERT_EQUALS(svi+svi,temp);
	}
	// Note: Multiplying two sparse vectors is undefined by design.
	void testMultiplication(void) {
		TS_ASSERT_EQUALS((2*svc)[0],4);
		TS_ASSERT_EQUALS((2*svs)[0],4);
		TS_ASSERT_EQUALS((2*svi)[0],4);
		TS_ASSERT_EQUALS((2*svl)[0],4);
		TS_ASSERT_EQUALS((2*svf)[0],4);
		TS_ASSERT_EQUALS((2*svd)[0],4);
	}
};

class gridTestSuite : public CxxTest::TestSuite
{
private:

public:
	void setUp(int argc, char* argv[]) {
		MMSP::Init(argc, argv);
	}
	void tearDown() {
		MMSP::Finalize();
	}
	void testSize1D() {
		MMSP::grid<1,MMSP::scalar<char> >   g1c(0, 0, 64);
		MMSP::grid<1,MMSP::scalar<short> >  g1s(0, 0, 64);
		MMSP::grid<1,MMSP::scalar<int> >    g1i(0, 0, 64);
		MMSP::grid<1,MMSP::scalar<long> >   g1l(0, 0, 64);
		MMSP::grid<1,MMSP::scalar<float> >  g1f(0, 0, 64);
		MMSP::grid<1,MMSP::scalar<double> > g1d(0, 0, 64);

		#ifndef MPI_VERSION
		TS_ASSERT_EQUALS(MMSP::nodes(g1c),64);
		TS_ASSERT_EQUALS(MMSP::nodes(g1s),64);
		TS_ASSERT_EQUALS(MMSP::nodes(g1i),64);
		TS_ASSERT_EQUALS(MMSP::nodes(g1l),64);
		TS_ASSERT_EQUALS(MMSP::nodes(g1f),64);
		TS_ASSERT_EQUALS(MMSP::nodes(g1d),64);
		#else
		int np = MPI::COMM_WORLD.Get_size();
		TS_ASSERT_LESS_THAN_EQUALS(MMSP::nodes(g1c),64/np);
		TS_ASSERT_LESS_THAN_EQUALS(MMSP::nodes(g1s),64/np);
		TS_ASSERT_LESS_THAN_EQUALS(MMSP::nodes(g1i),64/np);
		TS_ASSERT_LESS_THAN_EQUALS(MMSP::nodes(g1l),64/np);
		TS_ASSERT_LESS_THAN_EQUALS(MMSP::nodes(g1f),64/np);
		TS_ASSERT_LESS_THAN_EQUALS(MMSP::nodes(g1d),64/np);
		#endif

		TS_ASSERT_EQUALS(g1c.buffer_size(),MMSP::nodes(g1c)*sizeof(char));
		TS_ASSERT_EQUALS(g1s.buffer_size(),MMSP::nodes(g1s)*sizeof(short));
		TS_ASSERT_EQUALS(g1i.buffer_size(),MMSP::nodes(g1i)*sizeof(int));
		TS_ASSERT_EQUALS(g1l.buffer_size(),MMSP::nodes(g1l)*sizeof(long));
		TS_ASSERT_EQUALS(g1f.buffer_size(),MMSP::nodes(g1f)*sizeof(float));
		TS_ASSERT_EQUALS(g1d.buffer_size(),MMSP::nodes(g1d)*sizeof(double));

	}
};
