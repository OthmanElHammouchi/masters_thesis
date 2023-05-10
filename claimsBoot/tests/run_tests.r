if (require("RUnit", quietly = TRUE)) {
    require("patternBreak")
    test.dir <- system.file("unit_tests", package = "patternBreak")
    test.triangle <- ChainLadder::UKMotor

    testsuite <- defineTestSuite("patternBreak",
        dirs = test.dir,
        testFileRegexp = "^runit.+\\.r",
        testFuncRegexp = "^test.+",
        rngKind = "Mersenne-Twister",
        rngNormalKind = "Inversion")

    testResult <- runTestSuite(testsuite, verbose = 1)
    printTextProtocol(testResult)
}