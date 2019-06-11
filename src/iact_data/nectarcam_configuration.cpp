/*

   calin/iact_data/nectarcam_configuration.cpp -- Stephen Fegan -- 2017-12-13

   NectarCAM specific configuration

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <stdexcept>
#include <memory>

#include <calin_global_config.hpp>
#include <iact_data/nectarcam_configuration.hpp>
#include <util/log.hpp>
#include <util/string.hpp>
#include <util/file.hpp>

using namespace calin::util::log;

#ifdef CALIN_HAVE_XERCESC

#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>

using namespace xercesc;

namespace {

std::string to_string(const XMLCh* const to_transcode) {
  char* x = XMLString::transcode(to_transcode);
  std::string s(x);
  XMLString::release(&x);
  return s;
}

#if 0
std::string to_lower_string(const XMLCh* const to_transcode) {
  char* x = XMLString::transcode(to_transcode);
  std::string s(x);
  XMLString::release(&x);
  return calin::util::string::to_lower(s);
}
#endif

} // anonymous namespace

calin::ix::iact_data::nectarcam_configuration::NectarCamCameraConfiguration*
calin::iact_data::nectarcam_configuration::decode_nmc_xml_file(
  const std::string& filename,
  calin::ix::iact_data::nectarcam_configuration::NectarCamCameraConfiguration* nccc)
{
  calin::ix::iact_data::nectarcam_configuration::NectarCamCameraConfiguration* nccc_out = nullptr;

  try {
    XMLPlatformUtils::Initialize();
  }
  catch (const XMLException& toCatch) {
    LOG(WARNING) << "decode_nmc_xml_file: XERCES exception: " << to_string(toCatch.getMessage());
    return nullptr;
  }

  static const XMLCh gLS[] = { chLatin_L, chLatin_S, chNull };
  DOMImplementation *impl = DOMImplementationRegistry::getDOMImplementation(gLS);
  DOMLSParser       *parser = ((DOMImplementationLS*)impl)->createLSParser(DOMImplementationLS::MODE_SYNCHRONOUS, 0);
  DOMConfiguration  *config = parser->getDomConfig();

  config->setParameter(XMLUni::fgDOMNamespaces, false);
  config->setParameter(XMLUni::fgXercesSchema, false);
  config->setParameter(XMLUni::fgXercesHandleMultipleImports, true);
  config->setParameter(XMLUni::fgXercesSchemaFullChecking, false);
  config->setParameter(XMLUni::fgDOMDisallowDoctype, true);
  config->setParameter(XMLUni::fgDOMValidate, false);
  config->setParameter(XMLUni::fgDOMDatatypeNormalization, true);

  // And create our error handler and install it
  //DOMCountErrorHandler errorHandler;
  //config->setParameter(XMLUni::fgDOMErrorHandler, &errorHandler);

  if(not calin::util::file::is_readable(filename)) {
    LOG(WARNING) << filename << ": file does not exist";
    XMLPlatformUtils::Terminate();
    return nullptr;
  }

  try {
    DOMDocument* doc = parser->parseURI(filename.c_str());
    if (doc) {
      auto* doc_element = doc->getDocumentElement();
      while(doc_element and to_string(doc_element->getTagName())!="nectar") {
        doc_element = doc_element->getNextElementSibling();
      }
      if(doc_element) {
        nccc_out = nccc;
        if(nccc_out == nullptr)nccc_out = new calin::ix::iact_data::nectarcam_configuration::NectarCamCameraConfiguration;

        nccc_out->set_nmc_source_filename(filename);

        auto* element = doc_element->getFirstElementChild();
        while(element) {
          if(to_string(element->getTagName()) == "drawer") {
            auto* module = nccc_out->add_module();

            static const XMLCh att_number[] = { chLatin_n, chLatin_u, chLatin_m, chLatin_b, chLatin_e, chLatin_r, chNull };
            if(element->hasAttribute(att_number)) {
              module->set_module_id(calin::util::string::int32_from_string(to_string(element->getAttribute(att_number))));
            } else {
              module->set_module_id(-1);
            }

            static const XMLCh att_ip[] = { chLatin_i, chLatin_p, chNull };
            if(element->hasAttribute(att_ip)) {
              module->set_ip_address(to_string(element->getAttribute(att_ip)));
            }

            static const XMLCh att_mac[] = { chLatin_m, chLatin_a, chLatin_c, chNull };
            if(element->hasAttribute(att_mac)) {
              module->set_mac_address(to_string(element->getAttribute(att_mac)));
            }

            auto* sub_element = element->getFirstElementChild();
            while(sub_element) {
              if(to_string(sub_element->getTagName()) == "HVPA") {
                static const XMLCh att_voltage[7][9] = {
                  { chLatin_v, chLatin_o, chLatin_l, chLatin_t, chLatin_a, chLatin_g, chLatin_e, chDigit_0, chNull },
                  { chLatin_v, chLatin_o, chLatin_l, chLatin_t, chLatin_a, chLatin_g, chLatin_e, chDigit_1, chNull },
                  { chLatin_v, chLatin_o, chLatin_l, chLatin_t, chLatin_a, chLatin_g, chLatin_e, chDigit_2, chNull },
                  { chLatin_v, chLatin_o, chLatin_l, chLatin_t, chLatin_a, chLatin_g, chLatin_e, chDigit_3, chNull },
                  { chLatin_v, chLatin_o, chLatin_l, chLatin_t, chLatin_a, chLatin_g, chLatin_e, chDigit_4, chNull },
                  { chLatin_v, chLatin_o, chLatin_l, chLatin_t, chLatin_a, chLatin_g, chLatin_e, chDigit_5, chNull },
                  { chLatin_v, chLatin_o, chLatin_l, chLatin_t, chLatin_a, chLatin_g, chLatin_e, chDigit_6, chNull }
                };
                for(unsigned i=0; i<7; i++) {
                  if(sub_element->hasAttribute(att_voltage[i])) {
                    module->add_hvpa_voltage(calin::util::string::double_from_string(to_string(sub_element->getAttribute(att_voltage[i]))));
                  } else {
                    module->add_hvpa_voltage(-1.0);
                  }
                }
              } else if (to_string(sub_element->getTagName()) == "CntrlDaq") {
                static const XMLCh att_nf[] = { chLatin_n, chLatin_f, chNull };
                if(sub_element->hasAttribute(att_nf)) {
                  module->set_num_samples(calin::util::string::int32_from_string(to_string(sub_element->getAttribute(att_nf))));
                }

                static const XMLCh att_daqType[] = { chLatin_d, chLatin_a, chLatin_q, chLatin_T, chLatin_y, chLatin_p, chLatin_e, chNull };
                if(sub_element->hasAttribute(att_daqType)) {
                  module->set_daq_mode(to_string(sub_element->getAttribute(att_daqType)));
                }
              } else if (to_string(sub_element->getTagName()) == "CntrlFpga") {
                static const XMLCh att_NectarFreq[] = {
                  chLatin_N, chLatin_e, chLatin_c, chLatin_t, chLatin_a, chLatin_r,
                  chLatin_F, chLatin_r, chLatin_e, chLatin_q, chNull };
                if(sub_element->hasAttribute(att_NectarFreq)) {
                  uint32_t divisor = 1+calin::util::string::int32_from_string(to_string(sub_element->getAttribute(att_NectarFreq)));
                  module->set_clock_divisor(divisor);
                  module->set_nominal_sampling_frequency(2000.0/divisor);
                }
              }
              sub_element = sub_element->getNextElementSibling();
            }
          }
          element = element->getNextElementSibling();
        }
      }
    } else {
      LOG(WARNING) << filename << ": no document element found";
    }

    parser->release();
  }
  catch (const XMLException& toCatch) {
    LOG(WARNING) << filename << ": XERCES exception: " << to_string(toCatch.getMessage());
    parser->release();
    XMLPlatformUtils::Terminate();
  }
  catch (const DOMException& toCatch) {
    LOG(WARNING) << filename << ": parsing error: " << to_string(toCatch.getMessage());
    parser->release();
    XMLPlatformUtils::Terminate();
  }
  catch (...) {
    LOG(WARNING) << filename << ": unknown exception";
    parser->release();
    XMLPlatformUtils::Terminate();
  }

  return nccc_out;
}

#else

calin::ix::iact_data::nectarcam_configuration::NectarCamCameraConfiguration*
calin::iact_data::nectarcam_configuration::decode_nmc_xml_file(
  const std::string& filename,
  calin::ix::iact_data::nectarcam_configuration::NectarCamCameraConfiguration* nccc)
{
  LOG(WARNING) << "decode_nmc_xml_file: not implemented as XERCES not present at compile time.";
  return nullptr;
}

#endif
