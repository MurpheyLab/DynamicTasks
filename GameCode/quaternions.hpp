
double QuaternionLookup(int euler_angle, int pos_val) {
  double return_val=0.0;
  //printf("euler angle: %d.\n", euler_angle);
  if (euler_angle==90) {
      double quat_list[4] = {0, 0, 0.7071068, 0.7071068};
      return_val = quat_list[pos_val];
  } else if (euler_angle==100) {
      double quat_list[4] = {0, 0, 0.7660444, 0.6427876};
      return_val = quat_list[pos_val];
  } else if (euler_angle==110) {
      double quat_list[4] = {0, 0, 0.819152, 0.5735764};
      return_val = quat_list[pos_val];
  } else if (euler_angle==120) {
      double quat_list[4] = {0, 0, 0.8660254, 0.5};
      return_val = quat_list[pos_val];
  } else if (euler_angle==130) {
      double quat_list[4] = {0, 0, 0.9063078, 0.4226183};
      return_val = quat_list[pos_val];
  } else if (euler_angle==140) {
      double quat_list[4] = {0, 0, 0.9396926, 0.3420201};
      return_val = quat_list[pos_val];
  } else if (euler_angle==150) {
      double quat_list[4] = {0, 0, 0.9659258, 0.258819};
      return_val = quat_list[pos_val];
  } else if (euler_angle==160) {
      double quat_list[4] = {0, 0, 0.9848078, 0.1736482};
      return_val = quat_list[pos_val];
  } else if (euler_angle==170) {
      double quat_list[4] = {0, 0, 0.9961947, 0.0871557};
      return_val = quat_list[pos_val];
  } else if (euler_angle==180) {
      double quat_list[4] = {0, 0, 1, 0};
      return_val = quat_list[pos_val];
  } else if (euler_angle==190) {
      double quat_list[4] = {0, 0, 0.9961947, -0.0871557};
      return_val = quat_list[pos_val];
  } else if (euler_angle==200) {
      double quat_list[4] = {0, 0, 0.9848078, -0.1736482};
      return_val = quat_list[pos_val];
  } else if (euler_angle==210) {
      double quat_list[4] = {0, 0, 0.9659256230107514, -0.2588198037463362};
      return_val = quat_list[pos_val];
  } else if (euler_angle==220) {
      double quat_list[4] = {0, 0, 0.9396926, -0.3420201};
      return_val = quat_list[pos_val];
  } else if (euler_angle==230) {
      double quat_list[4] = {0, 0, 0.9063078, -0.4226183};
      return_val = quat_list[pos_val];
  } else if (euler_angle==240) {
      double quat_list[4] = {0, 0, 0.8660254, -0.5};
      return_val = quat_list[pos_val];
  } else if (euler_angle==250) {
      double quat_list[4] = {0, 0, 0.819152, -0.5735764};
      return_val = quat_list[pos_val];
  } else if (euler_angle==260) {
      double quat_list[4] = {0, 0, 0.7660444, -0.6427876};
      return_val = quat_list[pos_val];
  } else if (euler_angle==270) {
      double quat_list[4] = {0, 0, 0.7071068, -0.7071068};
      return_val = quat_list[pos_val];
  } else if (euler_angle == 225) {
      double quat_list[4] = { 0, 0, 0.92387921048041, -0.3826842098154747 };
      return_val = quat_list[pos_val];
  }
  return (return_val);
}
