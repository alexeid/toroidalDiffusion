// This is an empty umbrella build including all the component builds.
// This build is not necessarily needed. The component builds work independently.

rootProject.name = "toroidalDiffusion"

include("lphy")
include("lphy-studio")

pluginManagement {
    // the repos to load Gradle plugins
    repositories {
        mavenCentral()
//        maven {
//            // to local build/plugins
//            url = uri("${rootDir.parent}/GradlePlugins/platforms/build/releases/")
//            println("Temp repo : ${url}")
//        }
        // add sonatype snapshots repository
        maven {
            url=uri("https://s01.oss.sonatype.org/content/repositories/snapshots/")
        }
        gradlePluginPortal()
    }
}

