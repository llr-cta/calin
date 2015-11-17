/* 

   calin/unit_tests/io/test_vs_optics.cpp -- Stephen Fegan -- 2015-09-16

   Unit tests for VSOptics

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

   This file is part of "calin"
   
   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.
    
   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <gtest/gtest.h>
#include <google/protobuf/text_format.h>

#include <io/sql_transceiver.hpp>
#include <io/sqlite3_transceiver.hpp>
#include <simulation/vso_obscuration.hpp>
#include <simulation/vso_pixel.hpp>
#include <simulation/vso_mirror.hpp>
#include <simulation/vso_telescope.hpp>

using namespace calin::simulation::vs_optics;
using namespace calin::ix::simulation::vs_optics;
using namespace calin::io::sql_transceiver;

static double sqrt_eps = std::sqrt(std::numeric_limits<double>::epsilon());

TEST(TestVSOObscuration, StoreAndRetreive_Tube) {
  // Create tube obscuration and dump to proto
  VSOTubeObscuration obs1({1,2,3},{4,5,6},7,true);
  VSOObscurationData obs1_data;
  obs1.dump_to_proto(&obs1_data);

  // Store proto in SQL DB
  SQLite3Transceiver xvr("test_db.sqlite", SQLite3Transceiver::TRUNCATE_RW);
  //SQLite3Transceiver xvr(":memory:", SQLite3Transceiver::TRUNCATE_RW);
  xvr.create_tables("obscuration", VSOObscurationData::descriptor(),
                    nullptr, "Instance of obscuration");
  uint64_t oid;
  xvr.insert("obscuration", oid, &obs1_data);

  // Retreive proto from SQL DB
  VSOObscurationData obs2_data;
  xvr.retrieve_by_oid("obscuration", oid, &obs2_data);

  // Create obscuration from proto
  VSOObscuration* obs2_base = VSOObscuration::create_from_proto(obs2_data);
  VSOTubeObscuration* obs2 = dynamic_cast<VSOTubeObscuration*>(obs2_base);
  
  ASSERT_NE(obs2, nullptr);
  EXPECT_EQ(obs1.end1_pos(), obs2->end1_pos());
  EXPECT_EQ(obs1.end2_pos(), obs2->end2_pos());
  EXPECT_EQ(obs1.diameter(), obs2->diameter());
  EXPECT_EQ(obs1.incoming_only(), obs2->incoming_only());

  delete obs2;
}

TEST(TestVSOObscuration, StoreAndRetreive_Disk) {
  // Create tube obscuration and dump to proto
  VSODiskObscuration obs1({11,12,13},{14,15,16},17,true);
  VSOObscurationData obs1_data;
  obs1.dump_to_proto(&obs1_data);

  // Store proto in SQL DB
  SQLite3Transceiver xvr("test_db.sqlite", SQLite3Transceiver::EXISTING_RW);
  //SQLite3Transceiver xvr(":memory:", SQLite3Transceiver::TRUNCATE_RW);
  //xvr.create_tables("obscuration", VSOObscurationData::descriptor(),
  //nullptr, "Instance of obscuration");
  uint64_t oid;
  xvr.insert("obscuration", oid, &obs1_data);

  // Retreive proto from SQL DB
  VSOObscurationData obs2_data;
  xvr.retrieve_by_oid("obscuration", oid, &obs2_data);

  // Create obscuration from proto
  VSOObscuration* obs2_base = VSOObscuration::create_from_proto(obs2_data);
  VSODiskObscuration* obs2 = dynamic_cast<VSODiskObscuration*>(obs2_base);
  
  ASSERT_NE(obs2, nullptr);
  EXPECT_EQ(obs1.center_pos(), obs2->center_pos());
  EXPECT_EQ(obs1.normal(), obs2->normal());
  EXPECT_EQ(obs1.diameter(), obs2->diameter());
  EXPECT_EQ(obs1.incoming_only(), obs2->incoming_only());

  delete obs2;
}

TEST(TestVSOPixel, StoreAndRetreive) {
  // Create pixel and dump to proto
  VSOPixel pix1(nullptr,201,202,false,{203,204,205});
  VSOPixelData pix1_data;
  pix1.dump_to_proto(&pix1_data);

  // Store proto in SQL DB
  SQLite3Transceiver xvr("test_db.sqlite", SQLite3Transceiver::EXISTING_RW);
  //SQLite3Transceiver xvr(":memory:", SQLite3Transceiver::TRUNCATE_RW);
  xvr.create_tables("pixel", VSOPixelData::descriptor(),
                    nullptr, "Instance of pixel");
  uint64_t oid;
  xvr.insert("pixel", oid, &pix1_data);

  // Retreive proto from SQL DB
  VSOPixelData pix2_data;
  xvr.retrieve_by_oid("pixel", oid, &pix2_data);

  // Create pixel from proto
  VSOPixel* pix2 = VSOPixel::create_from_proto(pix2_data, nullptr);

  EXPECT_EQ(pix1.pos(), pix2->pos());
  EXPECT_EQ(pix1.removed(), pix2->removed());
  EXPECT_EQ(pix1.id(), pix2->id());
  EXPECT_EQ(pix1.hexID(), pix2->hexID());

  delete pix2;
}

TEST(TestVSOMirror, StoreAndRetreive) {
  VSOMirror mir1(nullptr,301,302,true,{303,304,305},{306,307,308},309,310,311);
  VSOMirrorData mir1_data;
  mir1.dump_to_proto(&mir1_data);

  // Store proto in SQL DB
  SQLite3Transceiver xvr("test_db.sqlite", SQLite3Transceiver::EXISTING_RW);
  //SQLite3Transceiver xvr(":memory:", SQLite3Transceiver::TRUNCATE_RW);
  xvr.create_tables("mirror", VSOMirrorData::descriptor(),
                    nullptr, "Instance of mirror");
  uint64_t oid;
  xvr.insert("mirror", oid, &mir1_data);

  // Retreive proto from SQL DB
  VSOMirrorData mir2_data;
  xvr.retrieve_by_oid("mirror", oid, &mir2_data);

  // Create mirror from proto
  VSOMirror* mir2 = VSOMirror::create_from_proto(mir2_data, nullptr);

  EXPECT_EQ(mir1.pos(), mir2->pos());
  EXPECT_EQ(mir1.align(), mir2->align());
  EXPECT_EQ(mir1.removed(), mir2->removed());
  EXPECT_EQ(mir1.id(), mir2->id());
  EXPECT_EQ(mir1.hexID(), mir2->hexID());
  EXPECT_EQ(mir1.focalLength(), mir2->focalLength());
  EXPECT_EQ(mir1.spotSize(), mir2->spotSize());
  EXPECT_EQ(mir1.degradingFactor(), mir2->degradingFactor());

  delete mir2;
}

TEST(TestVSOTelescope, StoreAndRetreive) {
  VSOTelescope scope1(1,2,{3,4,5},6,7,8,9,10,{11,12,13},14,15,16,17,18,19,
                      20,21,{22,23,24},25,26,27,28,29,{30,31,32},33,34);
#if 0
  for(unsigned i=0;i<100; i++)
    scope1.add_mirror(new VSOMirror(&scope1,100000.0+i,101000.0+i,true,
                                    {103000.0+i,104000.0+i,105000.0+i},
                                    {106000.0+i,107000.0+i,108000.0+i},
                                    109000.0+i,110000.0+i,111000.0+i));
  for(unsigned i=0;i<100; i++)
    scope1.add_pixel(new VSOPixel(&scope1,200000.0+i,201000.0+i,true,
                                  {203000.0+i,204000.0+i,205000.0+i}));
#endif
  for(unsigned i=0;i<50; i++)
    scope1.add_obscuration(
        new VSOTubeObscuration({300000.0+i,301000.0+i,302000.0+i},
                               {303000.0+i,304000.0+i,305000.0+i},
                               307000+i,true));
  for(unsigned i=50;i<100; i++)
    scope1.add_obscuration(
        new VSODiskObscuration({300000.0+i,301000.0+i,302000.0+i},
                               {303000.0+i,304000.0+i,305000.0+i},
                               307000+i,true));
  
  VSOTelescopeData scope1_data;
  scope1.dump_to_proto(&scope1_data);

  // Store proto in SQL DB
  SQLite3Transceiver xvr("test_db.sqlite", SQLite3Transceiver::EXISTING_RW,true);
  //SQLite3Transceiver xvr(":memory:", SQLite3Transceiver::TRUNCATE_RW);
  xvr.create_tables("scope", VSOTelescopeData::descriptor(),
                    nullptr, "Instance of scope");
  uint64_t oid;
  xvr.insert("scope", oid, &scope1_data);

  // Retreive proto from SQL DB
  VSOTelescopeData scope2_data;
  xvr.retrieve_by_oid("scope", oid, &scope2_data);

  for(unsigned i=0;i<scope1_data.obscuration_size();i++)
    std::cout << scope1_data.obscuration(i).has_disk() << ' '
              << scope1_data.obscuration(i).has_tube() << '\n';

  for(unsigned i=0;i<scope2_data.obscuration_size();i++)
    std::cout << scope2_data.obscuration(i).has_disk() << ' '
              << scope2_data.obscuration(i).has_tube() << '\n';

  // Create scope from proto
  VSOTelescope* scope2 = VSOTelescope::create_from_proto(scope2_data);

  EXPECT_EQ(scope1.id(), scope2->id());
  EXPECT_EQ(scope1.hexID(), scope2->hexID());
  EXPECT_EQ(scope1.altitude(), scope2->altitude());
  EXPECT_EQ(scope1.pos(), scope2->pos());
  EXPECT_NEAR(scope1.deltaY(), scope2->deltaY(), sqrt_eps);
  EXPECT_NEAR(scope1.alphaX(), scope2->alphaX(), sqrt_eps);
  EXPECT_NEAR(scope1.alphaY(), scope2->alphaY(), sqrt_eps);
  EXPECT_NEAR(scope1.azimuth(), scope2->azimuth(), sqrt_eps);
  EXPECT_NEAR(scope1.elevation(), scope2->elevation(), sqrt_eps);
  EXPECT_NEAR(scope1.opticalAxis().x, scope2->opticalAxis().x, sqrt_eps);
  EXPECT_NEAR(scope1.opticalAxis().y, scope2->opticalAxis().y, sqrt_eps);
  EXPECT_NEAR(scope1.opticalAxis().z, scope2->opticalAxis().z, sqrt_eps);

  EXPECT_EQ(scope1.translation(), scope2->translation());
  EXPECT_EQ(scope1.curvatureRadius(), scope2->curvatureRadius());
  EXPECT_EQ(scope1.aperture(), scope2->aperture());
  EXPECT_EQ(scope1.facetSpacing(), scope2->facetSpacing());
  EXPECT_EQ(scope1.facetSize(), scope2->facetSize());
  EXPECT_NEAR(scope1.reflectorRotation(), scope2->reflectorRotation(), sqrt_eps);
  EXPECT_EQ(scope1.mirrorHexRings(), scope2->mirrorHexRings());
  EXPECT_EQ(scope1.reflectorIP(), scope2->reflectorIP());
  EXPECT_EQ(scope1.mirrorParity(), scope2->mirrorParity());
  
  EXPECT_EQ(scope1.focalPlanePosition(), scope2->focalPlanePosition());
  EXPECT_EQ(scope1.cameraDiameter(), scope2->cameraDiameter());
  EXPECT_EQ(scope1.fov(), scope2->fov());
  EXPECT_EQ(scope1.cathodeDiameter(), scope2->cathodeDiameter());
  EXPECT_EQ(scope1.pixelSpacing(), scope2->pixelSpacing());
  EXPECT_EQ(scope1.concentratorSurvivalProb(),
            scope2->concentratorSurvivalProb());
  EXPECT_NEAR(scope1.focalPlaneRotion().x,scope2->focalPlaneRotion().x,sqrt_eps);
  EXPECT_NEAR(scope1.focalPlaneRotion().y,scope2->focalPlaneRotion().y,sqrt_eps);
  EXPECT_NEAR(scope1.focalPlaneRotion().z,scope2->focalPlaneRotion().z,sqrt_eps);
  EXPECT_EQ(scope1.cameraIP(), scope2->cameraIP());
  EXPECT_EQ(scope1.pixelParity(), scope2->pixelParity());

  EXPECT_EQ(scope1.numObscurations(), scope2->numObscurations());

  EXPECT_EQ(scope1.numMirrors(), scope2->numMirrors());
  EXPECT_EQ(scope1.numMirrorHexSites(), scope2->numMirrorHexSites());

  EXPECT_EQ(scope1.numPixels(), scope2->numPixels());
  EXPECT_EQ(scope1.numPixelHexSites(), scope2->numPixelHexSites());

  delete scope2;
}
 
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
