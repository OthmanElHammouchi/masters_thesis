if (require("RUnit", quietly = TRUE)) {

    require("patternBreak")

    test.dir <- system.file("unit_tests", package = "patternBreak")

    testsuite <- defineTestSuite("patternBreak",
        dirs = test.dir,
        testFileRegexp = "^runit.+\\.r",
        testFuncRegexp = "^test.+",
        rngKind = "Mersenne-Twister",
        rngNormalKind = "Inversion")

    testResult <- runTestSuite(testsuite, verbos = 0)
    printTextProtocol(testResult)

}