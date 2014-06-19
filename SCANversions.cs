﻿/* 
 * [Scientific Committee on Advanced Navigation]
 * 			S.C.A.N. Satellite
 *
 * SCANversions - logs version numbers for SCANsat and various associated mods
 * 
 * Copyright (c)2014 David Grandy <david.grandy@gmail.com>;
 * Copyright (c)2014 technogeeky <technogeeky@gmail.com>;
 * Copyright (c)2014 (Your Name Here) <your email here>; see LICENSE.txt for licensing details.
 */


using System.Collections.Generic;
using System.Reflection;
using System.Diagnostics;
using System.Linq;
using UnityEngine;

namespace SCANsat
{
    [KSPAddon(KSPAddon.Startup.MainMenu, true)]
    internal class SCANversions: MonoBehaviour
    {
        private string[] Assemblies = new string[7] {"SCANsat", "SCANsatRPM", "SCANsatKethane", "Kethane", "RasterPropMonitor", "MechJebRPM", "MechJeb2"};
        internal static string SCANsatVersion = "";
        internal static string SCANurl = "";
        private List<AssemblyLog> assemblyList = new List<AssemblyLog>();

        private void Start() {
            findAssemblies(Assemblies);
        }

        private void findAssemblies(string[] assemblies) {
            foreach (string name in assemblies) { //Search for the relevant plugins among the loaded assemblies
                AssemblyLoader.LoadedAssembly assembly = null;
                try {
                    assembly = AssemblyLoader.loadedAssemblies.SingleOrDefault(a => a.assembly.GetName().Name == name);
                    if (assembly != null)
                        assemblyList.Add(new AssemblyLog(assembly));
                }
                catch {
                    assembly = AssemblyLoader.loadedAssemblies.FirstOrDefault(a => a.assembly.GetName().Name == name);
                    if (assembly != null) {
                        assemblyList.Add(new AssemblyLog(assembly));
                        ScreenMessages.PostScreenMessage(string.Format("Multiple copies of assembly: {0}.dll detected, only one copy should be installed", name), 15f, ScreenMessageStyle.UPPER_CENTER);
                    }
                }
            }
            if (assemblyList.Count > 0) { 
                SCANsatVersion = assemblyList[0].infoVersion;
                SCANurl = assemblyList[0].location;
                debugWriter();
            }
        }

        private void debugWriter() {
            foreach (AssemblyLog log in assemblyList) {
                print(string.Format("[SCANlogger] Assembly: {0} found; Version: {1}; Informational Version: {2}; Location: {3}", log.name, log.version, log.infoVersion, log.location));
            }           
        }

    }

    //A class to gather and store information about assemblies
    internal class AssemblyLog
    {
        internal string name, version, infoVersion, location;

        internal AssemblyLog(AssemblyLoader.LoadedAssembly assembly)
        {
            name = assembly.assembly.GetName().Name;
            version = assembly.assembly.GetName().Version.ToString();
            infoVersion = FileVersionInfo.GetVersionInfo(assembly.assembly.Location).ProductVersion; 
            location = assembly.url.ToString();
        }
    
    }


}
