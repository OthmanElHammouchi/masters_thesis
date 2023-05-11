if (require("RUnit", quietly = TRUE)) {
    require("claimsBoot")
    test.dir <- system.file("unit_tests", package = "claimsBoot")
    test.triangle <- ChainLadder::UKMotor

    testsuite <- defineTestSuite("claimsBoot",
        dirs = test.dir,
        testFileRegexp = "^runit.+\\.r",
        testFuncRegexp = "^test.+",
        rngKind = "Mersenne-Twister",
        rngNormalKind = "Inversion")

    testResult <- runTestSuite(testsuite, verbose = 1)
    printTextProtocol(testResult)
}