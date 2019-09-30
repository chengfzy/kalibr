// THIS FILE HAS BEEN AUTOGENERATED BY gen_files.py
#include <aslam/Frame.hpp>
#include <aslam/FrameBaseSerialization.hpp>
#include <aslam/cameras.hpp>

// Standard serialization headers
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
// These ones are in sm_boost
#include <boost/portable_binary_iarchive.hpp>
#include <boost/portable_binary_oarchive.hpp>

BOOST_CLASS_EXPORT_IMPLEMENT(aslam::Frame<aslam::cameras::MaskedOmniCameraGeometry>);

namespace aslam {

template void Frame<aslam::cameras::MaskedOmniCameraGeometry>::save<>(boost::archive::text_oarchive& ar,
                                                                      const unsigned int version) const;
template void Frame<aslam::cameras::MaskedOmniCameraGeometry>::load<>(boost::archive::text_iarchive& ar,
                                                                      const unsigned int version);
template void Frame<aslam::cameras::MaskedOmniCameraGeometry>::save<>(boost::archive::xml_oarchive& ar,
                                                                      const unsigned int version) const;
template void Frame<aslam::cameras::MaskedOmniCameraGeometry>::load<>(boost::archive::xml_iarchive& ar,
                                                                      const unsigned int version);
template void Frame<aslam::cameras::MaskedOmniCameraGeometry>::save<>(boost::archive::binary_oarchive& ar,
                                                                      const unsigned int version) const;
template void Frame<aslam::cameras::MaskedOmniCameraGeometry>::load<>(boost::archive::binary_iarchive& ar,
                                                                      const unsigned int version);
template void Frame<aslam::cameras::MaskedOmniCameraGeometry>::save<>(boost::archive::portable_binary_oarchive& ar,
                                                                      const unsigned int version) const;
template void Frame<aslam::cameras::MaskedOmniCameraGeometry>::load<>(boost::archive::portable_binary_iarchive& ar,
                                                                      const unsigned int version);

}  // namespace aslam
