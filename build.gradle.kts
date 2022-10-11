// buildSrc is a trap, use composite builds
// https://docs.gradle.org/current/userguide/structuring_software_products.html

// this is the umbrella build to define cross-build lifecycle tasks.
// https://docs.gradle.org/current/userguide/structuring_software_products_details.html

import java.text.SimpleDateFormat
import java.util.*

plugins {
    `java-library`
//    `maven-publish`
}

// Configures this project and each of its sub-projects.
allprojects {
    repositories {
        mavenCentral()
        // add sonatype snapshots repository
        maven {
            url=uri("https://s01.oss.sonatype.org/content/repositories/snapshots/")
        }
        // Managing plugin versions via pluginManagement in settings.gradle.kts
//        mavenLocal() // only for testing
    }

}

// disable build folder at root project
tasks.jar {enabled = false}
tasks.build {enabled = false}

// Configures the sub-projects of this project.
subprojects {
    group = "io.github.linguaphylo"
    version = "0.0.1"//-SNAPSHOT"
    val webSteam = "github.com/alexeid/toroidalDiffusion"
    val web = "https://${webSteam}"
    val homepage = web // may be different to github page

    var calendar: Calendar? = Calendar.getInstance()
    var formatter = SimpleDateFormat("dd-MMM-yyyy HH:mm:ss")

    // shared attributes
    tasks.withType<Jar>() {
        manifest {
            attributes(
                "Implementation-Version" to archiveVersion,
                "Implementation-URL" to web,
                "Implementation-Vendor" to "Alexei Drummond",
                "Built-By" to "Alexei Drummond", //System.getProperty("user.name"),
                "Build-Jdk" to JavaVersion.current().majorVersion.toInt(),
                "Built-Date" to formatter.format(calendar?.time)
            )
        }
        // copy LICENSE to META-INF
        metaInf {
            from(rootDir) {
                include("LICENSE")
            }
        }

    }

}

