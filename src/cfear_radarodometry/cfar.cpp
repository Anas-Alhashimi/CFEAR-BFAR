#include "cfear_radarodometry/cfar.h"

// BFAR //////////////////////// old implementation
void BFAR_filter(cv_bridge::CvImagePtr &polar, pcl::PointCloud<pcl::PointXYZI>::Ptr& cloud, int window_size_, double scale_factor, double offset_factor_, double range_res, double min_distance)
{
	if(cloud==NULL)
		cloud = pcl::PointCloud<pcl::PointXYZI>::Ptr (new pcl::PointCloud<pcl::PointXYZI>());

	float scan_thr_mat[400][4000] = {{0}};
	float scan_float[400][4000] = {{0}};
	bool  det_mat[400][4000] = {{0}};
    // added by Anas to convert the scan into float and to calculate azimuth mean
    // float azimuth_mean[400] = {0};
	float scan_mean = 0;
    for (int bearing = 0; bearing < polar->image.rows; bearing++){
		for (size_t i = 0; i < polar->image.cols; i++){
			scan_float[bearing][i] = (float)(polar->image.at<uchar>(bearing, i));
			scan_float[bearing][i] = pow(10.0f,scan_float[bearing][i]/40.0f);
			//azimuth_mean[bearing]  = azimuth_mean[bearing] + scan_float[bearing][i];
			scan_mean = scan_mean + scan_float[bearing][i];
		}
		//azimuth_mean[bearing] = azimuth_mean[bearing] /(polar->image.cols);
	}
	scan_mean = scan_mean/(float)(polar->image.rows*polar->image.cols);
	
	
	int nb_guard_cells_ = 2;
	
	float trailing_mean, mean, forwarding_mean, max_value_thr;
	int trailing_window_start, trailing_window_end, forwarding_window_start, forwarding_window_end;
	for (int azimuth_nb = 0; azimuth_nb < polar->image.rows; azimuth_nb++){
		for(int range_bin = 0; range_bin < polar->image.cols; range_bin++){
			float cells_number = 0.0;
			trailing_window_start = std::max(0, range_bin - nb_guard_cells_ - window_size_);
			trailing_window_end = range_bin - nb_guard_cells_;
			//trailing_mean = getMean(azimuth, trailing_window_start, trailing_window_end);
			float sum = 0.;
			//int N = 0;
			for(int i = trailing_window_start; i < trailing_window_end; i++){
				sum += scan_float[azimuth_nb][i];
				cells_number += 1.;
				}
			trailing_mean = sum;
			//cells_number = (float)N;
			forwarding_window_start = range_bin + nb_guard_cells_;
			forwarding_window_end = std::min(polar->image.cols, range_bin + nb_guard_cells_ + window_size_);
			//forwarding_mean = getMean(azimuth, forwarding_window_start, forwarding_window_end);
			sum = 0.;
			//N = 0;
			for(int i = forwarding_window_start; i < forwarding_window_end; i++){
				sum += scan_float[azimuth_nb][i];
				cells_number += 1.;
				}
			forwarding_mean = sum;
			//cells_number = cells_number + (float)N;
			mean = (trailing_mean + forwarding_mean)/cells_number;//(2.0 * (float)window_size_);
			//cout << cells_number <<" " ;
			
			scan_thr_mat[azimuth_nb][range_bin] = offset_factor_ + mean * scale_factor;
			//cout <<"offset_factor_"<<"="<<offset_factor_<<" --";
			//std::cout <<"scale_factor"<<"="<<scale_factor<<" --";
			det_mat[azimuth_nb][range_bin] = (scan_float[azimuth_nb][range_bin] > scan_thr_mat[azimuth_nb][range_bin]);
		  }
		}
  pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_nofilter(new pcl::PointCloud<pcl::PointXYZI>);
  const  double min_distance_sqrd = min_distance*min_distance;
  sensor_msgs::ImagePtr msg = polar->toImageMsg();
  float theta;
  /*if(cv_polar_image->image.rows!=400 || cv_polar_image->image.cols!=3768){
    std::cout<<"Size error rows: "<<cv_polar_image->image.rows<<", cols:"<<cv_polar_image->image.cols<<std::endl;
    exit(0);
  }*/
  for (int bearing = 0; bearing < polar->image.rows; bearing++){
    theta = ((float)(bearing+1) / polar->image.rows) * 2 * M_PI;
    std::vector<pcl::PointXYZI> pnts_sorted;
    for (size_t i = 0; i < polar->image.cols; i++){
      int ind = i+0; //Unsure about this one!
      //double d = range_res*(ind);

      /*if(d < min_distance || d > max_distance){
        continue;
      }*/
      pcl::PointXYZI p;
      p.x = range_res * ind * cos(theta);
      p.y = range_res * ind * sin(theta);
      p.intensity = polar->image.at<uchar>(bearing, ind);
      p.z = 0.0;//(p.intensity-z_min)/10.0;
      
      if (det_mat[bearing][ind] == 1 ) // only push the good ones :)
		pnts_sorted.push_back(p);
      cloud_nofilter->push_back(p);
      //InsertStrongestK(pnts_sorted, p, k_strongest);
    }
    for(auto && p :pnts_sorted){
      if (p.x*p.x+p.y*p.y>min_distance_sqrd)
        cloud->push_back(p);
    }
  }

  cloud->width = (int)cloud->points.size();
  cloud->height = 1;



  /*cloud_nofilter->width = (int)cloud_nofilter->points.size();
  cloud_nofilter->height = 1;
  cloud_nofilter->header.frame_id = radar_frameid;*/

  pcl_conversions::toPCL(polar->header.stamp, cloud->header.stamp);//pcl_conversions::toPCL(cv_polar_image->header.stamp,cloud->header.stamp);
  //pcl_conversions::toPCL(cv_polar_image->header.stamp, cloud_nofilter->header.stamp);//pcl_conversions::toPCL(cv_polar_image->header.stamp,cloud->header.stamp);
  //FilteredPublisher.publish(cloud);
  //UnfilteredPublisher.publish(cloud_nofilter);
}
// end of BFAR

/////////////////////////////////// CFAR Filter ////////////////////////////////////////

CFARFilter::CFARFilter(const double &false_alarm_rate,  const double &range_resolution,
                       const double &static_threshold, const double &min_distance, const double &max_distance) :
  false_alarm_rate_(false_alarm_rate), range_resolution_(range_resolution), static_threshold_(static_threshold),
  min_distance_(min_distance), max_distance_(max_distance)
{}

double CFARFilter::getCAScalingFactor(const int &window_size) const
{
  const double N = window_size;
  return  N * (std::pow(false_alarm_rate_, -1./N) - 1.);
}

double CFARFilter::getIntensity(uchar intensity, bool square_law) const
{
  if(square_law)
    return std::pow(double(intensity), 2.);

  return double(intensity);
}

//////////////////////////////////// Azimuth CA-CFAR ////////////////////////////////////

AzimuthCACFAR::AzimuthCACFAR(const int &window_size, const double &false_alarm_rate, const int &nb_guard_cells,
                             const double &range_resolution, const double &static_threshold, const double &min_distance, const double &max_distance) :
  CFARFilter(false_alarm_rate, range_resolution, static_threshold, min_distance, max_distance), window_size_(window_size), nb_guard_cells_(nb_guard_cells)
{
  this->scaling_factor_ = getCAScalingFactor(window_size_ * 2);
}

void AzimuthCACFAR::getFilteredPointCloud(const cv_bridge::CvImagePtr &radar_image, pcl::PointCloud<pcl::PointXYZI>::Ptr &output_pointcloud) const
{
  for(int azimuth_nb = 0; azimuth_nb < radar_image->image.rows; azimuth_nb++)
  {
    cv::Mat azimuth = radar_image->image.row(azimuth_nb);
    const double theta = (double(azimuth_nb + 1) / radar_image->image.rows) * 2. * M_PI;
    for(int range_bin = 0; range_bin < azimuth.cols; range_bin++)
    {
      const double range = range_resolution_ * double(range_bin);
      const double intensity = double(azimuth.at<uchar>(range_bin));
      if(range > min_distance_ && range < max_distance_ && intensity > static_threshold_) //static th not officially part of CA-CFAR but speeds up and makes result more accurate. This
      {
        //Scaling factor should be updated with actual window size.
        const int trailing_window_start = std::max(0, range_bin - nb_guard_cells_ - window_size_);
        const int trailing_window_end = range_bin - nb_guard_cells_;
        const double trailing_mean = getMean(azimuth, trailing_window_start, trailing_window_end);

        const int forwarding_window_start = range_bin + nb_guard_cells_;
        const int forwarding_window_end = std::min(azimuth.cols, range_bin + nb_guard_cells_ + window_size_);
        const double forwarding_mean = getMean(azimuth, forwarding_window_start, forwarding_window_end);

        const double mean = (trailing_mean + forwarding_mean)/2.0;    // CA-CFAR
        //const double mean std::max(trailing_mean, forwarding_mean); // GO-CFAR
        const double threshold = scaling_factor_ * mean;
        const double squared_intensity = std::pow(intensity, 2.);
        if(squared_intensity > threshold)
        {
          pcl::PointXYZI p;
          p.x = range * std::cos(theta);
          p.y = range * std::sin(theta);
          p.intensity = intensity;
          output_pointcloud->push_back(p);
        }
      }
    }
  }
}

double AzimuthCACFAR::getMean(const cv::Mat &azimuth, const int &start_idx, const int &end_idx) const
{
  double sum = 0.;
  double N = 0.;
  for(size_t i = start_idx; i < end_idx; i++)
  {
    sum += std::pow(double(azimuth.at<uchar>(i)), 2.);
    N += 1.;
  }
  return sum / N;
}

//////////////////////////////////// Azimuth CA-BFAR //////////////////////////////////// Anas

AzimuthBFAR::AzimuthBFAR(const int &window_size, const double &scale_factor, const double &offset_factor, const int &nb_guard_cells,
                             const double &range_resolution, const double &min_distance, const double &max_distance) : 
  CFARFilter(1, range_resolution, 1, min_distance, max_distance), window_size_(window_size), nb_guard_cells_(nb_guard_cells),offset_factor_(offset_factor),scale_factor_(scale_factor)
{
  //this->scaling_factor_ = getCAScalingFactor(window_size_ * 2);
}
void AzimuthBFAR::getFilteredPointCloud(const cv_bridge::CvImagePtr &radar_image, pcl::PointCloud<pcl::PointXYZI>::Ptr &output_pointcloud) const
{
  for(int azimuth_nb = 0; azimuth_nb < radar_image->image.rows; azimuth_nb++)
  {
    cv::Mat azimuth = radar_image->image.row(azimuth_nb);
    const double theta = (double(azimuth_nb + 1) / radar_image->image.rows) * 2. * M_PI;
    for(int range_bin = 0; range_bin < azimuth.cols; range_bin++)
    {
      const double range = range_resolution_ * double(range_bin);
      const double intensity = double(azimuth.at<uchar>(range_bin));
      if(range > min_distance_ && range < max_distance_ ) // Anas static th not officially part of CA-CFAR but speeds up and makes result more accurate. This
      {
        
        const int trailing_window_start = std::max(0, range_bin - nb_guard_cells_ - window_size_);
        const int trailing_window_end = range_bin - nb_guard_cells_;
        const double trailing_mean = getMean(azimuth, trailing_window_start, trailing_window_end);
		//std::cout<< window_size_<<",";
        const int forwarding_window_start = range_bin + nb_guard_cells_;
        const int forwarding_window_end = std::min(azimuth.cols, range_bin + nb_guard_cells_ + window_size_);
        const double forwarding_mean = getMean(azimuth, forwarding_window_start, forwarding_window_end);

        const double mean = (trailing_mean + forwarding_mean)/2.0;    // CA-CFAR
        //const double mean std::max(trailing_mean, forwarding_mean); // GO-CFAR
        const double threshold = offset_factor_ + scale_factor_ * mean;
        const double scaled_intensity = pow(10.0f,intensity/40.0f);
        if(scaled_intensity > threshold)
        {
          pcl::PointXYZI p;
          p.x = range * std::cos(theta);
          p.y = range * std::sin(theta);
          p.intensity = intensity;
          output_pointcloud->push_back(p);
        }
      }
    }
  }
}

double AzimuthBFAR::getMean(const cv::Mat &azimuth, const int &start_idx, const int &end_idx) const
{
  double sum = 0.;
  double N = 0.;
  for(size_t i = start_idx; i < end_idx; i++)
  {
    sum += std::pow(double(azimuth.at<uchar>(i)), 2.);
    N += 1.;
  }
  return sum / N;
}
