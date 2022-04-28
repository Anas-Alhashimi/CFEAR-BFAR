#include "cfear_radarodometry/odometrykeyframefuser.h"
namespace CFEAR_Radarodometry {

visualization_msgs::Marker GetDefault(){
  visualization_msgs::Marker m;
  m.color.r = 0;
  m.color.g = 0;
  m.color.b = 1;
  m.color.a = 1;
  m.type = visualization_msgs::Marker::LINE_LIST;
  m.id = 0;
  m.scale.x = 0.1;
  m.action = visualization_msgs::Marker::ADD;
  m.lifetime = ros::Duration(0);
  m.header.frame_id ="world";
  m.header.stamp = ros::Time::now();
  return m;
}

OdometryKeyframeFuser::OdometryKeyframeFuser(const Parameters& pars, bool disable_callback) : par(pars), nh_("~"){
  assert (!par.input_points_topic.empty() && !par.scan_registered_latest_topic.empty() && !par.scan_registered_keyframe_topic.empty() && !par.odom_latest_topic.empty() && !par.odom_keyframe_topic.empty() );
  assert(par.res>0.05 && par.submap_scan_size>=1 );
  //radar_reg =
  radar_reg = boost::shared_ptr<CFEAR_Radarodometry::n_scan_normal_reg>(new n_scan_normal_reg(Str2Cost(par.cost_type),
                                                                                              Str2loss(par.loss_type_),
                                                                                              par.loss_limit_,
                                                                                              par.weight_opt));

  radar_reg->SetD2dPar(par.covar_scale_,par.regularization_);

  Tprev_fused = Eigen::Affine3d::Identity();
  Tcurrent = Eigen::Affine3d::Identity();
  T_prev = Eigen::Affine3d::Identity();
  Tmot = Eigen::Affine3d::Identity();

  pose_current_publisher = nh_.advertise<nav_msgs::Odometry>(par.odom_latest_topic,50);
  pose_keyframe_publisher = nh_.advertise<nav_msgs::Odometry>(par.odom_keyframe_topic,50);
  pubsrc_cloud_latest = nh_.advertise<pcl::PointCloud<pcl::PointXYZI> >(par.scan_registered_latest_topic, 1000);
  pub_cloud_keyframe = nh_.advertise<pcl::PointCloud<pcl::PointXYZI> >(par.scan_registered_keyframe_topic, 1000);

  if(!disable_callback) {
    cout<<"subscribe<sensor_msgs::PointCloud2>("<<par.input_points_topic<<")"<<endl;
    pointcloud_callback = nh_.subscribe<sensor_msgs::PointCloud2>(par.input_points_topic, 1000,
                                                                  &OdometryKeyframeFuser::pointcloudCallback, this,
                                                                  ros::TransportHints().tcpNoDelay(true));
  }
  else
    cout<<"callback disabled"<<endl;
}
/*OdometryKeyframeFuser::~OdometryKeyframeFuser(){
  //cout<<"destruct fuser"<<endl;
  keyframes_.clear();
  ///  cout<<"destructed fuser"<<endl;
}*/



bool OdometryKeyframeFuser::KeyFrameBasedFuse(const Eigen::Affine3d& diff, bool use_keyframe, double min_keyframe_dist, double min_keyframe_rot_deg){

  if(!use_keyframe)
    return true;
  Eigen::Vector3d Tmotion_euler = diff.rotation().eulerAngles(0,1,2);
  CFEAR_Radarodometry::normalizeEulerAngles(Tmotion_euler);
  bool fuse_frame = false;
  //cout<<"diff (trans[m]/rot[deg])=("<<diff.translation().norm()<<"/"<<diff.rotation().eulerAngles(0,1,2).norm()*180.0/M_PI<<") limit=("<<min_keyframe_dist_<<"/"<<min_keyframe_rot_deg_<<")"<<endl;
  if(diff.translation().norm()>min_keyframe_dist || Tmotion_euler.norm()>(min_keyframe_rot_deg*M_PI/180.0))
    fuse_frame = true;
  return fuse_frame;
}


bool OdometryKeyframeFuser::AccelerationVelocitySanityCheck(const Eigen::Affine3d& Tmot_prev, const Eigen::Affine3d& Tmot_curr) {
  const double dt = 0.25 ; // 4Hz
  const double vel_limit = 200; // m/s
  const double acc_limit = 200; // m/s
  const double vel = (Tmot_curr.translation()/dt).norm();
  const double acc = ((Tmot_curr.translation() - Tmot_prev.translation())/(dt*dt)).norm();

  //cout<<"vel: "<<vel<<", acc: "<<acc<<endl;
  if(acc>acc_limit){
    cout<<"acceleration exceeds limit"<<acc<<" > "<<acc_limit<<endl;
    return false;
  }
  else if(vel>vel_limit){
    cout<<"velocity exceeds limit"<<vel<<" > "<<vel_limit<<endl;
    return false;
  }
  else return true;

}


Eigen::Affine3d OdometryKeyframeFuser::Interpolate(const Eigen::Affine3d &T2, double factor, const Eigen::Affine3d &T1) {
  // std::cout<<"factor="<<factor<<std::endl;
  Eigen::Affine3d pose;
  Eigen::Quaterniond q1(T1.linear());
  Eigen::Quaterniond q2(T2.linear());
  pose.linear() = (q1.slerp(factor, q2)).toRotationMatrix();
  Eigen::Vector3d tdiff = T2.translation() - T1.translation();
  pose.translation() = T1.translation() + tdiff * factor;
  return pose;
}
pcl::PointXYZI OdometryKeyframeFuser::Transform(const Eigen::Affine3d& T, pcl::PointXYZI& p){
  Eigen::Vector3d v(p.x,p.y,p.z);
  Eigen::Vector3d v2 = T*v;
  pcl::PointXYZI p2;
  p2.x = v2(0);
  p2.y = v2(1);
  p2.z = v2(2);
  p2.intensity = p.intensity;
  return p2;
}
nav_msgs::Odometry OdometryKeyframeFuser::FormatOdomMsg(const Eigen::Affine3d& T, const Eigen::Affine3d& Tmot, const ros::Time& t, Matrix6d& Cov){
  nav_msgs::Odometry odom_msg;

  //double d = Tmot.translation().norm();
  Eigen::MatrixXd cov_1_36(Cov);
  //reg_cov(0,0) = reg_cov(1,1) = (d*0.02)*(d*0.02);
  //reg_cov(5,5) = (d*0.005)*(d*0.005);
  cov_1_36.resize(1,36);
  for(int i=0;i<36;i++){
    odom_msg.pose.covariance[i] = cov_1_36.data()[i];
  }

  odom_msg.header.stamp = t;
  odom_msg.header.frame_id = par.odometry_link_id;
  odom_msg.child_frame_id = "sensor";
  tf::poseEigenToMsg( T, odom_msg.pose.pose);
  return odom_msg;
}
pcl::PointCloud<pcl::PointXYZI> OdometryKeyframeFuser::FormatScanMsg(pcl::PointCloud<pcl::PointXYZI>& cloud_in, Eigen::Affine3d& T){
  pcl::PointCloud<pcl::PointXYZI> cloud_out;
  pcl::transformPointCloud(cloud_in, cloud_out, T);
  cloud_out.header.frame_id = par.odometry_link_id;
  cloud_out.header.stamp = cloud_in.header.stamp;
  return cloud_out;
}

void OdometryKeyframeFuser::processFrame(pcl::PointCloud<pcl::PointXYZI>::Ptr& cloud) {

  ros::Time t0 = ros::Time::now();
  if(par.compensate)
    Compensate(*cloud, Tmot, par.radar_ccw);

  ros::Time t;
  pcl_conversions::fromPCL(cloud->header.stamp, t);

  std::vector<Matrix6d> cov_vek;
  std::vector<CFEAR_Radarodometry::MapNormalPtr> scans_vek;
  std::vector<Eigen::Affine3d> T_vek;
  ros::Time t1 = ros::Time::now();
  CFEAR_Radarodometry::MapNormalPtr Pcurrent = CFEAR_Radarodometry::MapNormalPtr(new MapPointNormal(cloud, par.res, Eigen::Vector2d(0,0), par.weight_intensity_, par.use_raw_pointcloud));





  ros::Time t2 = ros::Time::now();
  Eigen::Affine3d Tguess;
  if(par.use_guess)
    Tguess = T_prev*Tmot;
  else
    Tguess = T_prev;

  if(keyframes_.empty()){
    AddToReference(keyframes_, Pcurrent, Eigen::Affine3d::Identity(), 2);
    return;
  }
  else
    FormatScans(keyframes_, Pcurrent, Tguess, cov_vek, scans_vek, T_vek);


  bool success = true;
  Eigen::Affine3d Tprior = T_prev;
  //cout<<"prior:"<<Tprior.translation().transpose()<<endl;
  //cout<<"guess:"<<Tguess.translation().transpose()<<endl;

  if(!par.disable_registration)
    bool success = radar_reg->Register(scans_vek, T_vek, cov_vek, par.soft_constraint);
    //bool success = radar_reg->RegisterTimeContinuous(scans_vek, T_vek, cov_vek, Tmot, par.soft_constraint, par.radar_ccw);



    //bool success = radar_reg->RegisterTimeContinuous(scans_vek, T_vek, cov_vek, Tmot, par.soft_constraint, par.radar_ccw);
    //bool success = radar_reg->Register(scans_vek, T_vek, cov_vek, par.soft_constraint);
    //

  ros::Time t3 = ros::Time::now();
  //cout<<radar_reg->getScore()<<endl;

  if(success==false){
    //    cout<<"registration failure"<<radar_reg->summary_.FullReport()<<endl;
    exit(0);
  }

  //PlotAssociations(scans_vek, T_vek, radar_reg->scan_associations_,  t);
  Tcurrent = T_vek.back();
  //Pcurrent->Compensate(Tmot, par.radar_ccw);
  Eigen::Affine3d Tmot_current = T_prev.inverse()*Tcurrent;
  if(!AccelerationVelocitySanityCheck(Tmot, Tmot_current))
    Tcurrent = Tguess; //Insane acceleration and speed, lets use the guess.
  Tmot = T_prev.inverse()*Tcurrent;


  MapPointNormal::PublishMap("/current_normals", Pcurrent, Tcurrent, par.odometry_link_id,-1);
  pcl::PointCloud<pcl::PointXYZI> cld_latest = FormatScanMsg(*cloud, Tcurrent);
  nav_msgs::Odometry msg_current = FormatOdomMsg(Tcurrent, Tmot, t, cov_vek.back());
  pubsrc_cloud_latest.publish(cld_latest);
  pose_current_publisher.publish(msg_current);
  if(par.publish_tf_){
    geometry_msgs::TransformStamped transformStamped;
    transformStamped.header.stamp = msg_current.header.stamp;
    transformStamped.header.frame_id = msg_current.header.frame_id;

    tf::Transform Tf;
    std::vector<tf::StampedTransform> trans_vek;
    tf::transformEigenToTF(Tcurrent, Tf);
    trans_vek.push_back(tf::StampedTransform(Tf, t, par.odometry_link_id, "sensor_est"));
    trans_vek.push_back(tf::StampedTransform(Tf, t, par.odometry_link_id, "base_link"));
    Tbr.sendTransform(trans_vek);
  }

  const Eigen::Affine3d Tkeydiff = keyframes_.back().first.inverse()*Tcurrent;
  bool fuse = KeyFrameBasedFuse(Tkeydiff, par.use_keyframe, par.min_keyframe_dist_, par.min_keyframe_rot_deg_);

  CFEAR_Radarodometry::timing.Document("velocity", Tmot.translation().norm()/Tsensor);

  /*ndt_map::NDTMapMsg msg_ndt;
  //static void ToNDTMsg(std::vector<cellptr>& cells, ndt_map::NDTMapMsg& msg);
  //auto cells = Pcurrent->TransformCells(Tcurrent);

  auto cells = Pcurrent->GetCells();
  cell::ToNDTMsg(cells,msg_ndt);
  msg_ndt.header.stamp = t;
  msg_ndt.header.frame_id = "sensor_est";
  static ros::Publisher pub_map = nh_.advertise<ndt_map::NDTMapMsg>("/ndt_map",100);
  pub_map.publish(msg_ndt);*/

  if(success && fuse){
    distance_traveled += Tkeydiff.translation().norm();
    Tprev_fused = Tcurrent;
    pcl::PointCloud<pcl::PointXYZI> cld_keyframe = FormatScanMsg(*cloud, Tcurrent);
    nav_msgs::Odometry msg_keyframe = FormatOdomMsg(Tcurrent, Tkeydiff, t, cov_vek.back());
    pub_cloud_keyframe.publish(cld_keyframe);
    pose_keyframe_publisher.publish(msg_keyframe);

    frame_nr_++;
    AddToReference(keyframes_, Pcurrent, Tprev_fused, 2);
  }
  ros::Time t4 = ros::Time::now();
  CFEAR_Radarodometry::timing.Document("compensate", ToMs(t1-t0));
  CFEAR_Radarodometry::timing.Document("build_normals", ToMs(t2-t1));
  CFEAR_Radarodometry::timing.Document("register", ToMs(t3-t2));
  CFEAR_Radarodometry::timing.Document("publish_etc", ToMs(t4-t3));
  T_prev = Tcurrent;
}

void OdometryKeyframeFuser::pointcloudCallback(const sensor_msgs::PointCloud2::ConstPtr& msg_in){

  //  cout<<"callback"<<endl;
  ros::Time t = ros::Time::now();
  pcl::PointCloud<pcl::PointXYZI>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZI>());
  pcl::fromROSMsg (*msg_in, *cloud);
  pcl_conversions::toPCL(msg_in->header.stamp, cloud->header.stamp);
  this->processFrame(cloud);
  nr_callbacks_++;
  ros::Time t2 = ros::Time::now();
  CFEAR_Radarodometry::timing.Document("Registration-full",CFEAR_Radarodometry::ToMs(t2-t));
}


void OdometryKeyframeFuser::pointcloudCallback(const pcl::PointCloud<pcl::PointXYZI>::Ptr& msg_in, Eigen::Affine3d &Tcurr){
  ros::Time t = ros::Time::now();
  pcl::PointCloud<pcl::PointXYZI>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZI>());
  *cloud = *msg_in;
  cloud->header = msg_in->header;
  this->processFrame(cloud);
  nr_callbacks_++;
  Tcurr = Tcurrent;
  ros::Time t2 = ros::Time::now();
  CFEAR_Radarodometry::timing.Document("Registration",CFEAR_Radarodometry::ToMs(t2-t));

}

void OdometryKeyframeFuser::PrintSurface(const std::string& path, const Eigen::MatrixXd& surface){
  std::ofstream myfile;
  myfile.open (path);
  for(int i=0;i<surface.rows();i++){
    for(int j=0;j<surface.cols();j++){
      if(j==surface.cols()-1)
        myfile<<surface(i,j)<<std::endl;
      else
        myfile<<surface(i,j)<<" ";
    }
  }
  myfile.close();
}

void AddToReference(PoseScanVector& reference, MapNormalPtr cloud,  const Eigen::Affine3d& T, size_t submap_scan_size){
  reference.push_back( std::make_pair(T, cloud) );
  if(reference.size() > submap_scan_size){
    reference.erase(reference.begin());
  }
}

void FormatScans(const PoseScanVector& reference,
                                        const CFEAR_Radarodometry::MapNormalPtr& Pcurrent,
                                        const Eigen::Affine3d& Tcurrent,
                                        std::vector<Matrix6d>& cov_vek,
                                        std::vector<MapNormalPtr>& scans_vek,
                                        std::vector<Eigen::Affine3d>& T_vek
                                        ){

  for (int i=0;i<reference.size();i++) {
    cov_vek.push_back(Identity66);
    scans_vek.push_back(reference[i].second);
    T_vek.push_back(reference[i].first);
  }
  cov_vek.push_back(Identity66);
  scans_vek.push_back(Pcurrent);
  T_vek.push_back(Tcurrent);
}

}
