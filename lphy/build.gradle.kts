plugins {
    `java-library`
}

base.archivesName.set("toroidaldiffusion-lphy")

java {
    sourceCompatibility = JavaVersion.VERSION_17
    targetCompatibility = JavaVersion.VERSION_17
    withSourcesJar()
}

dependencies {
    api("io.github.linguaphylo:lphy:1.3.1")

    implementation("org.ejml:ejml-core:0.38")
    implementation("org.ejml:ejml-ddense:0.38")
    implementation("org.ejml:ejml-simple:0.38")
}

tasks.jar {
    manifest {
        // shared attr in the root build
        attributes(
            "Implementation-Title" to "Toroidal diffusion",
            "Implementation-Vendor" to "LPhy developer team",
        )
    }
}

tasks.test {
    useJUnit()
    // useJUnitPlatform()
    // set heap size for the test JVM(s)
    minHeapSize = "128m"
    maxHeapSize = "1G"
    // show standard out and standard error of the test JVM(s) on the console
    testLogging.showStandardStreams = true
    //testLogging.exceptionFormat = org.gradle.api.tasks.testing.logging.TestExceptionFormat.FULL
}

