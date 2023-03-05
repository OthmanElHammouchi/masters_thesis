if (require("RUnit", quietly = TRUE)) {

    require("patternBreak")

    load(file.path(system.file(package = "patternBreak"), "data", "test_triangle.RData"))

    test.dir <- system.file("unit_tests", package = "patternBreak")

    testsuite <- defineTestSuite("patternBreak",
        dirs = test.dir,
        testFileRegexp = "^runit.+\\.r",
        testFuncRegexp = "^test.+",
        rngKind = "Mersenne-Twister",
        rngNormalKind = "Inversion")

    testResult <- runTestSuite(testsuite, verbose = 1)
    printTextProtocol(testResult)

}