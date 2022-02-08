#pragma once

#include "pcl/point_cloud.h"
#include "pcl/point_types.h"
#include "pcl_ros/point_cloud.h"
#include "sensor_msgs/Image.h"
#include <cv_bridge/cv_bridge.h>
#include <image_transport/image_transport.h>
#include "ros/subscriber.h"
#include "ros/publisher.h"

#include "pcl_conversions/pcl_conversions.h"
#include "cfear_radarodometry/radar_filters.h"
#include "cfear_radarodometry/statistics.h"

namespace CFEAR_Radarodometry  {
using std::cout;
using std::endl;
using std::cerr;


class radarDriver
{
public:
  class Parameters
  {
  public:
    Parameters() {}

    float z_min = 60; // min power
    float range_res = 0.0438;
    int azimuths = 400, k_strongest = 12;
    float min_distance = 2.5, max_distance = 200;
    std::string radar_frameid = "sensor_est", topic_filtered = "/Navtech/Filtered";
    std::string dataset = "oxford";

    void GetParametersFromRos( ros::NodeHandle& param_nh){
      param_nh.param<float>("range_res", range_res, 0.0438);
      param_nh.param<float>("z_min", z_min, 60);
      param_nh.param<float>("min_distance", min_distance, 2.5);
      param_nh.param<float>("max_distance", max_distance, 130);
      param_nh.param<int>("kstrongest", k_strongest, 12);
      param_nh.param<std::string>("topic_filtered", topic_filtered, "/Navtech/Filtered");
      param_nh.param<std::string>("radar_frameid", radar_frameid, "sensor_est");
      param_nh.param<std::string>("dataset", dataset, "oxford");
    }
    std::string ToString(){
      std::ostringstream stringStream;
      //stringStream << "RadarDriver::Parameters:"<<endl;
      stringStream << "range res, "<<range_res<<endl;
      stringStream << "z min, "<<z_min<<endl;
      stringStream << "min distance, "<<min_distance<<endl;
      stringStream << "max distance, "<<max_distance<<endl;
      stringStream << "k strongest, "<<k_strongest<<endl;
      stringStream << "topic_filtered, "<<topic_filtered<<endl;
      stringStream << "radar_frameid, "<<radar_frameid<<endl;
      stringStream << "dataset, "<<dataset<<endl;

      return stringStream.str();
    }

  };

  radarDriver(const Parameters& pars, bool disable_callback = false);

  ~radarDriver(){cout<<"destruct driver"<<endl;}

  void CallbackOffline(const sensor_msgs::ImageConstPtr &radar_image_polar, pcl::PointCloud<pcl::PointXYZI>::Ptr &cloud);

private:

  void InitAngles();

  void Callback(const sensor_msgs::ImageConstPtr &radar_image_polar);

  void CallbackOxford(const sensor_msgs::ImageConstPtr &radar_image_polar);


  Parameters par;
  std::vector<float> sin_values;
  std::vector<float> cos_values;
  float max_distance_sqrd, min_distance_sqrd;
  ros::NodeHandle nh_;
  ros::Subscriber sub;
  ros::Publisher FilteredPublisher,ExperimentalPublisher, UnfilteredPublisher;
  image_transport::Publisher pub;
  image_transport::ImageTransport it;


  pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_filtered;


};

}
