/* MMSP.unitTest.hpp
 * Use Google Test to build unit tests for MMSP classes
 * http://github.com/google/googletest
 */

#include"MMSP.hpp"
#include<gtest/gtest.h>

// GTest supports templated unit testing. Declare the types of interest.

using testing::Types;
typedef Types<char,short,int,long,float,double> datatypes;

template <class T>
class dataTestSuite : public testing::Test
{
protected:
	virtual void SetUp() {
		myScalar = 2;

		MMSP::vector<T> temp(3,2);
		myVector = temp;

		for (int i=0; i<3; i++)
			mySparse.set(2*i) = 2;
	}
	MMSP::scalar<T> myScalar;
	MMSP::vector<T> myVector;
	MMSP::sparse<T> mySparse;
};

TYPED_TEST_CASE(dataTestSuite, datatypes);

TYPED_TEST(dataTestSuite, testSize) {
	/*scalar*/	EXPECT_EQ(this->myScalar.buffer_size(),sizeof(TypeParam));
	/*vector*/	EXPECT_EQ(this->myVector.length(),3);
	          	EXPECT_EQ(this->myVector.buffer_size(),sizeof(int)+3*sizeof(TypeParam));
	          	EXPECT_THROW(this->myVector[4],std::out_of_range);

	/*sparse*/	EXPECT_EQ(this->mySparse.length(),3);
	          	EXPECT_EQ(this->mySparse.buffer_size(),sizeof(int) + 6*std::max(sizeof(int),sizeof(TypeParam)));
}
TYPED_TEST(dataTestSuite, testValue) {
	/*scalar*/	EXPECT_EQ(this->myScalar,static_cast<TypeParam>(2));
	/*vector*/	for (int i=0; i<3; i++)
	          		EXPECT_EQ(this->myVector[i],static_cast<TypeParam>(2));
	/*sparse*/	for (int i=0; i<6; i+=2)
	          		EXPECT_EQ(this->mySparse[i],static_cast<TypeParam>(2));
	          	for (int i=1; i<7; i+=2)
	          		EXPECT_EQ(this->mySparse[i],0);
}
TYPED_TEST(dataTestSuite, testAddition) {
	/*scalar*/	EXPECT_EQ(this->myScalar+this->myScalar,static_cast<TypeParam>(4));
	/*vector*/	MMSP::vector<TypeParam> tempv = this->myVector+this->myVector;
	          	for (int i=0; i<3; i++)
	          		EXPECT_EQ(tempv[i],static_cast<TypeParam>(4));
	/*sparse*/	MMSP::sparse<TypeParam> temps = this->mySparse+this->mySparse;
	          	for (int i=0; i<3; i++)
	          		EXPECT_EQ(temps.value(i),static_cast<TypeParam>(4));
}
TYPED_TEST(dataTestSuite, testIncrement) {
	/*scalar*/	this->myScalar += this->myScalar;
	          	EXPECT_EQ(this->myScalar,static_cast<TypeParam>(4));
	/*vector*/	this->myVector += this->myVector;
	          	for (int i=0; i<3; i++)
	          		EXPECT_EQ(this->myVector[i],static_cast<TypeParam>(4));
	/*sparse*/	this->mySparse+=this->mySparse;
	          	for (int i=0; i<3; i++)
	          		EXPECT_EQ(this->mySparse.value(i),static_cast<TypeParam>(4));
}
TYPED_TEST(dataTestSuite, testSquares) {
	/*scalar*/	EXPECT_EQ(this->myScalar*this->myScalar,static_cast<TypeParam>(4));
	/*vector*/	EXPECT_EQ(this->myVector*this->myVector,12);
	/*sparse*/	// No operator*(sparse,sparse): meaning is not obvious.
}
TYPED_TEST(dataTestSuite, testMultiplication) {
	/*scalar*/	EXPECT_EQ(2*this->myScalar,static_cast<TypeParam>(4));
	          	EXPECT_EQ(this->myScalar*2,static_cast<TypeParam>(4));
	/*vector*/	MMSP::vector<TypeParam> tempv = 2*this->myVector;
	          	MMSP::vector<TypeParam> tempv2 = this->myVector*2;
	          	for (int i=0; i<3; i++) {
	          		EXPECT_EQ(tempv[i],static_cast<TypeParam>(4));
	          		EXPECT_EQ(tempv2[i],static_cast<TypeParam>(4));
	          	}
	/*sparse*/	MMSP::sparse<TypeParam> temps = 2*this->mySparse;
	          	MMSP::sparse<TypeParam> temps2 = this->mySparse*2;
	          	for (int i=0; i<3; i++) {
	          		EXPECT_EQ(temps.value(i),static_cast<TypeParam>(4));
	          		EXPECT_EQ(temps2.value(i),static_cast<TypeParam>(4));
	          	}
}


template<class T>
class gridTest : public testing::Test
{
protected:
	// MUST call grid constructor in test constructor
	gridTest() : grid1D(1, 0, 64) {
		for (int n=0; n<nodes(grid1D); n++)
			grid1D(n) = n;
	};

	virtual void SetUp(int argc, char* argv[]) {
		MMSP::Init(argc, argv);
	}
	void tearDown() {
		MMSP::Finalize();
	}

	MMSP::grid<1,MMSP::scalar<T> > grid1D;
};

TYPED_TEST_CASE(gridTest, datatypes);

TYPED_TEST(gridTest, testSize) {
	#ifndef MPI_VERSION
	EXPECT_EQ(MMSP::nodes(this->grid1D),64);
	#else
	int np = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	EXPECT_LE(MMSP::nodes(this->grid1D),1+64/np);
	#endif
	EXPECT_EQ(this->grid1D.buffer_size(),MMSP::nodes(this->grid1D)*sizeof(TypeParam));
}

