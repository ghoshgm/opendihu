#include "output_writer/exfile/exfile.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>
#include <ctime>
#include <iomanip>

#include <iostream>
#include <vector>

#include "easylogging++.h"

#include <utility/python_utility.h>
#include <utility/petsc_utility.h>
#include <mesh/structured_regular_fixed.h>
#include <mesh/structured_deformable.h>
#include <mesh/unstructured_deformable.h>
#include <mesh/mesh.h>

namespace OutputWriter
{

void Exfile::outputComFile()
{
  std::stringstream s;
  s <<filenameBase_<< ".com";
  std::string filenameCom = s.str();
  // open file
  std::ofstream file = openFile(filenameCom);

  // filename without path
  std::string basename = filenameBase_;

  if (basename.rfind("/") != std::string::npos)
  {
    basename = basename.substr(basename.rfind("/")+1);
  }

  time_t time = std::time(NULL);
  struct tm *t = localtime(&time);

  file << "# Cmgui visualization file, generated by opendihu at "
    << t->tm_mday << "/" << t->tm_mon+1 << "/" << t->tm_year+1900 << " "
    << std::setfill('0') << std::setw(2) << t->tm_hour << ":" << t->tm_min << ":" << t->tm_sec << std::endl
    << "# run by cmgui " << basename << ".com" << std::endl << std::endl;
  file << "# read files" << std::endl;

  std::string sphereSize = this->specificSettings_.getOptionString("sphereSize", "0.005*0.005*0.01");
  
  std::string lastMeshName;
  int lastDimensionality = 3;
  
  // loop over times of stored files
  for (std::map<double,std::vector<FilenameWithElementAndNodeCount>>::iterator iter = filenamesWithElementAndNodeCount_.begin(); iter != filenamesWithElementAndNodeCount_.end(); iter++)
  {
    double currentTime = iter->first;
    std::vector<FilenameWithElementAndNodeCount> &files = iter->second;
      
    node_no_t nodeOffset = 1;
    element_no_t elementOffset = 1;
    
    int timeNo = 1;
    
    // loop over stored files with their contained number of nodes and elements
    for (std::vector<FilenameWithElementAndNodeCount>::iterator iter2 = files.begin(); iter2 != files.end(); iter2++)
    {
      // filename without path
      std::string basename = iter2->filename;

      if (basename.rfind("/") != std::string::npos)
      {
        basename = basename.substr(basename.rfind("/")+1);
      }

      file << "$fname = \"" << basename << "\"" << std::endl
        << "gfx read nodes node_offset " << nodeOffset << " time " << currentTime << " $fname.\".exnode\"" << std::endl
        << "gfx read elements node_offset " << nodeOffset << " line_offset 1 face_offset 1 element_offset " << elementOffset << " $fname.\".exelem\"" << std::endl << std::endl;
        //<< "gfx modify g_element $group surface;" << std::endl
        //<< "$n+=1500;" << std::endl;
        
      // increase offsets by number of nodes and elements in the current files
      nodeOffset += iter2->nNodes;
      elementOffset += iter2->nElements;
      lastMeshName = iter2->meshName;
      lastDimensionality = iter2->dimensionality;
      timeNo++;
    }
  }
  file << std::endl
    << "# set the group name of the mesh which should be output below. This is by default the last loaded dataset." << std::endl
    << "$group = \"" << lastMeshName << "\"" << std::endl;
    
  if (lastDimensionality == 2 || lastDimensionality == 3)
  {
    file << std::endl
      << "##### 3D mesh #####" << std::endl
      << "# create materials" << std::endl
      << "gfx create material muscle_transparent normal_mode ambient 0.4 0.14 0.11 diffuse 0.5 0.12 0.1 emission 0 0 0 specular 0.3 0.5 0.5 alpha 0.25 shininess 0.2" << std::endl << std::endl
      << "# create faces" << std::endl
      << "gfx define face egroup $group" << std::endl << std::endl
      << "# add line representation" << std::endl
      << "gfx modify g_element \"/\" lines" << std::endl << std::endl
      << "# add surface representation with muscle material" << std::endl
      << "gfx modify g_element $group surfaces material muscle_transparent" << std::endl << std::endl;
  }
  else if (lastDimensionality == 1)
  {
    file << std::endl
      << "##### 1D mesh #####" << std::endl
      << "# add spheres representation" << std::endl
      << "gfx modify g_element $group points domain_mesh1d coordinate geometry glyph sphere size \"" << sphereSize << "\" select_on material default data solution" << std::endl << std::endl;
  }
  
  file << "# add axes" << std::endl
    << "gfx modify g_element \"/\" points domain_point glyph axes_xyz size \"50*50*50\" select_on material silver selected_material default" << std::endl << std::endl
    << "# open Graphics Window" << std::endl
    << "gfx cre win" << std::endl << std::endl
    << "# open Scene Editor" << std::endl
    << "gfx edit scene" << std::endl << std::endl;

  file.close();
  
  LOG(INFO) << "File \"" << filenameCom << "\" written.";
}

};
