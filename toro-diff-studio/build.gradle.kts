plugins {
    `java-library`
}

//base.archivesName.set("toroidaldiffusion-lphy-studio")

java {
    sourceCompatibility = JavaVersion.VERSION_17
    targetCompatibility = JavaVersion.VERSION_17
    withSourcesJar()
}

dependencies {
    api("io.github.linguaphylo:lphy-studio:1.3.3-SNAPSHOT")

    // it is toroidalDiffusion/toro-diff-lphy
    implementation(project(":toro-diff-lphy"))

}

tasks.jar {
    manifest {
        // shared attr in the root build
        attributes(
            "Implementation-Title" to "Toroidal diffusion GUI"
        )
    }
}

//tasks.test {
//    useJUnit()
//    // useJUnitPlatform()
//    // set heap size for the test JVM(s)
//    minHeapSize = "128m"
//    maxHeapSize = "1G"
//    // show standard out and standard error of the test JVM(s) on the console
//    testLogging.showStandardStreams = true
//    //testLogging.exceptionFormat = org.gradle.api.tasks.testing.logging.TestExceptionFormat.FULL
//}

