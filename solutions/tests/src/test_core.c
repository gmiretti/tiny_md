#include "unity.h"
#include "core.h"

// every test file requires this function
// setUp() is called before each test case function
void setUp(void) {
    printf("Setting up stuff before test\n");
}

// every test file requires this function
// tearDown() is called after each test case function
void tearDown(void) {
    printf("Cleaning up stuff after test\n");
}

void test_function_should_doBlahAndBlah(void) {
    //test stuff
    printf("Ignoring that we are testing that we should do BlahAndBlah\n");
    TEST_IGNORE_MESSAGE("This Test Was Ignored On Purpose");
}

void test_function_should_doAlsoDoBlah(void) {
    //more test stuff
    printf("Testing that we should do also Blah\n");
}

int main(void) {
    UNITY_BEGIN();
    // We should list all tests here
    RUN_TEST(test_function_should_doBlahAndBlah);
    RUN_TEST(test_function_should_doAlsoDoBlah);
    return UNITY_END();
}