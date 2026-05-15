#!/bin/bash
echo "PREBUILD SCRIPT"
echo "Check to make sure that the meson compiler is using the conda environment."
echo "a common error is the compiler falls back to an incompatible system compiler."
echo "environment: python: $(which python)"
echo "Environment: gfortran: $(which gfortran)"
echo "Environment: f2py: $(which f2py)"
echo "..."
echo ""



rm -rf builddir

# export FC=$(which gfortran) 
# this explicity set the fortran compiler, 
# I used this when I had fotran as an argument for the project in the meson.build file "project('tstrippy', ['fortran','c'],version: '0.0.1',license: 'MIT', )"
# having fortran set here prompts meson to search for the fortran compiler based on the shell variable
# so by setting it here, I can use the fortran compiler in the meson.build file
# however, this sucks and I want to make sure that it is handeled by python and conda.
# so now, the commands below in the meson file select the fortran compiler based on the conda environment

# Configure
meson setup builddir 
    # --native-file <(echo "[binaries]"; echo "fortran = 'gfortran'")


# Build
meson compile -C builddir
meson install -C builddir/

echo ""

echo ""
echo "_________________________________________.____________________________.___. "
echo "\__    ___/   _____/\__    ___/\______   \   \______   \______   \__  |   | "
echo "  |    |  \_____  \   |    |    |       _/   ||     ___/|     ___//   |   | "
echo "  |    |  /        \  |    |    |    |   \   ||    |    |    |    \____   | "
echo "  |____| /_______  /  |____|    |____|_  /___||____|    |____|    / ______| "
echo "                 \/                    \/                         \/        "
echo "                         .___                                               "
echo "  ____   ____   ____   __| _/______                                         "
echo " /    \_/ __ \_/ __ \ / __ |/  ___/                                         "
echo "|   |  \  ___/\  ___// /_/ |\___ \                                          "
echo "|___|  /\___  >\___  >____ /____  >                                         "
echo "     \/     \/     \/     \/    \/                                          "
echo " .----------------.  .----------------.  .----------------.  "
echo "| .--------------. || .--------------. || .--------------. | "
echo "| |  ____  ____  | || |     ____     | || | _____  _____ | | "
echo "| | |_  _||_  _| | || |   .'    '.   | || ||_   _||_   _|| | "
echo "| |   \ \  / /   | || |  /  .--.  \  | || |  | |    | |  | | "
echo "| |    \ \/ /    | || |  | |    | |  | || |  | '    ' |  | | "
echo "| |    _|  |_    | || |  \  '--'  /  | || |   \ '--' /   | | "
echo "| |   |______|   | || |   '.____.'   | || |    '.__.'    | | "
echo "| |              | || |              | || |              | | "
echo "| '--------------' || '--------------' || '--------------' | "
echo " '----------------'  '----------------'  '----------------'  "
echo ""
echo "  __                        .__                                             "
echo "_/  |_  ____     __________ |  |___  __ ____                                "
echo "\   __\/  _ \   /  ___/  _ \|  |\  \/ // __ \                               "
echo " |  | (  <_> )  \___ (  <_> )  |_\   /\  ___/                               "
echo " |__|  \____/  /____  >____/|____/\_/  \___  >                              "
echo "                    \/                     \/                               "
echo "  ___ ___                .__.__   __              /\                        "
echo " /   |   \_____    _____ |__|  |_/  |_  ____   ___)/  ______                "
echo "/    ~    \__  \  /     \|  |  |\   __\/  _ \ /    \ /  ___/                "
echo "\    Y    // __ \|  Y Y  \  |  |_|  | (  <_> )   |  \\___ \                 "
echo " \___|_  /(____  /__|_|  /__|____/__|  \____/|___|  /____  >                "
echo "       \/      \/      \/                         \/     \/                 "
echo "___________                    __  .__                                      "
echo "\_   _____/ ________ _______ _/  |_|__| ____   ____   ______                "
echo " |    __)_ / ____/  |  \__  \\   __\  |/  _ \ /    \ /  ___/                "
echo " |        < <_|  |  |  // __ \|  | |  (  <_> )   |  \\___ \                 "
echo "/_______  /\__   |____/(____  /__| |__|\____/|___|  /____  >                "
echo "        \/    |__|          \/                    \/     \/                 "

echo "BUILD SUCCESS"