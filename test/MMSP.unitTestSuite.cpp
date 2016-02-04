/* MMSP.unitTest.hpp
** Use Google Test to build unit tests for MMSP classes
** http://github.com/google/googletest
*/
#include"MMSP.hpp"
#include<gtest/gtest.h>

using testing::Types;
typedef Types<char,short,int,long,float,double> datatypes;



template <class T>
class scalarTest : public testing::Test
{
protected:
	virtual void SetUp() {
		value = 2;
	}

	MMSP::scalar<T> value;
};

TYPED_TEST_CASE(scalarTest, datatypes);

TYPED_TEST(scalarTest, testSize) {
	EXPECT_EQ(this->value.buffer_size(),sizeof(TypeParam));
}
TYPED_TEST(scalarTest, testAddition) {
	EXPECT_EQ(this->value+this->value,static_cast<TypeParam>(4));
}
TYPED_TEST(scalarTest, testSquares) {
	EXPECT_EQ(this->value*this->value,static_cast<TypeParam>(4));
}
TYPED_TEST(scalarTest, testMultiplication) {
	EXPECT_EQ(2*this->value,static_cast<TypeParam>(4));
	EXPECT_EQ(this->value*2,static_cast<TypeParam>(4));
}




template <class T>
class vectorTest : public testing::Test
{
protected:
	virtual void SetUp() {
		MMSP::vector<T> temp(3,2);
		value = temp;
	}

	MMSP::vector<T> value;
};

TYPED_TEST_CASE(vectorTest, datatypes);

TYPED_TEST(vectorTest, testValue) {
	for (int i=0; i<3; i++) {
		EXPECT_EQ(this->value[i],static_cast<TypeParam>(2));
	}
}
TYPED_TEST(vectorTest, testSize) {
	EXPECT_EQ(this->value.length(),3);
	EXPECT_EQ(this->value.buffer_size(),sizeof(TypeParam));
	EXPECT_DEATH(this->value[4],"index exceeds vector size");
}
TYPED_TEST(vectorTest, testAddition) {
	MMSP::vector<TypeParam> temp = this->value+this->value;
	for (int i=0; i<3; i++)
		EXPECT_EQ(temp[i],static_cast<TypeParam>(4));
}
TYPED_TEST(vectorTest, testIncrement) {
	this->value+=this->value;
	for (int i=0; i<3; i++)
		EXPECT_EQ(this->value[i],static_cast<TypeParam>(4));
}
TYPED_TEST(vectorTest, testSquares) {
	EXPECT_EQ(this->value*this->value,12);
}
TYPED_TEST(vectorTest, testMultiplication) {
	MMSP::vector<TypeParam> a = 2*this->value;
	//MMSP::vector<TypeParam> b = this->value*2;
	for (int i=0; i<3; i++) {
		EXPECT_EQ(a[i],static_cast<TypeParam>(4));
		//EXPECT_EQ(b[i],static_cast<TypeParam>(4)); // no such operator
	}
}




template<class T>
class sparseTest : public testing::Test
{
protected:
	virtual void SetUp() {
		for (int i=0; i<3; i++)
			value.set(2*i) = 2;
	}

	MMSP::sparse<T> value;
};


TYPED_TEST_CASE(sparseTest, datatypes);

TYPED_TEST(sparseTest, testValue) {
	for (int i=0; i<6; i+=2)
		EXPECT_EQ(this->value[i],static_cast<TypeParam>(2));
	for (int i=1; i<7; i+=2)
		EXPECT_EQ(this->value[i],0);
}
TYPED_TEST(sparseTest, testSize) {
	EXPECT_EQ(this->value.length(),3);
	unsigned int delta = std::abs(sizeof(int) - sizeof(TypeParam));
	EXPECT_EQ(this->value.buffer_size(),sizeof(int) + 2*3*sizeof(TypeParam) + 3*delta);
}
TYPED_TEST(sparseTest, testAddition) {
	MMSP::sparse<TypeParam> temp = this->value+this->value;
	for (int i=0; i<3; i++)
		EXPECT_EQ(temp.value(i),static_cast<TypeParam>(4));
}
TYPED_TEST(sparseTest, testIncrement) {
	this->value+=this->value;
	for (int i=0; i<3; i++)
		EXPECT_EQ(this->value.value(i),static_cast<TypeParam>(4));
}
TYPED_TEST(sparseTest, testMultiplication) {
	MMSP::sparse<TypeParam> a = 2*this->value;
	//MMSP::sparse<TypeParam> b = this->value*2;
	for (int i=0; i<3; i++) {
		EXPECT_EQ(a.value(i),static_cast<TypeParam>(4));
		//EXPECT_EQ(b.value(i),static_cast<TypeParam>(4)); // no such operator
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
	int np = MPI::COMM_WORLD.Get_size();
	EXPECT_LE(MMSP::nodes(this->grid1D),64/np);
	#endif
	EXPECT_EQ(this->grid1D.buffer_size(),MMSP::nodes(this->grid1D)*sizeof(TypeParam));
}


