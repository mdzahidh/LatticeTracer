#define CONV_PPIPED_OLD_INPLACE_RegionOne_Even(value) \
   pdata[0] = Fetch(startX, startY, startZ); \
   pdata[1] = Fetch(startX, startY-1, startZ-1); \
   pdata[2] = Fetch(startX-1, startY, startZ-1); \
   pdata[3] = Fetch(startX-1, startY-1, startZ+1); \
   pdata[4] = Fetch(startX-1, startY, startZ); \
   pdata[5] = Fetch(startX, startY-1, startZ); \
   pdata[6] = Fetch(startX, startY, startZ-2); \
   pdata[7] = Fetch(startX-1, startY-1, startZ-1); \
   pdata[8] = Fetch(startX, startY, startZ+1); \
   pdata[9] = Fetch(startX, startY+1, startZ); \
   pdata[10] = Fetch(startX, startY, startZ+2); \
   pdata[11] = Fetch(startX+1, startY+1, startZ+2); \
   pdata[12] = Fetch(startX, startY, startZ+3); \
   pdata[13] = Fetch(startX, startY+1, startZ+1); \
   pdata[14] = Fetch(startX-1, startY, startZ+1); \
   pdata[15] = Fetch(startX, startY+1, startZ+2); \
   pdata[16] = Fetch(startX, startY, startZ-1); \
   pdata[17] = Fetch(startX+1, startY+1, startZ); \
   pdata[18] = Fetch(startX+1, startY, startZ-2); \
   pdata[19] = Fetch(startX, startY+1, startZ-2); \
   pdata[20] = Fetch(startX, startY, startZ-3); \
   pdata[21] = Fetch(startX, startY+1, startZ-1); \
   pdata[22] = Fetch(startX+1, startY, startZ-1); \
   pdata[23] = Fetch(startX+1, startY+1, startZ-2); \
   pdata[24] = Fetch(startX+1, startY, startZ); \
   pdata[25] = Fetch(startX, startY-1, startZ+1); \
   pdata[26] = Fetch(startX+1, startY, startZ+1); \
   pdata[27] = Fetch(startX+1, startY-1, startZ-1); \
   pdata[28] = Fetch(startX+2, startY, startZ); \
   pdata[29] = Fetch(startX+1, startY-1, startZ); \
   pdata[30] = Fetch(startX+1, startY, startZ+2); \
   pdata[31] = Fetch(startX+1, startY-1, startZ+1); \
   alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f; \
   value += (-10.0f*pdata[0]+4.0f*(pdata[1]+pdata[2]+pdata[3])-2.0f*(pdata[4]+pdata[5]+pdata[6])+pdata[7]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[0]-2.0f*(pdata[2]+pdata[3])+pdata[4]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[3])+pdata[5]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[2])+pdata[6]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[0]+pdata[3]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[2]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[1]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[0]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f; \
   value += (-10.0f*pdata[8]+4.0f*(pdata[9]+pdata[10]+pdata[11])-2.0f*(pdata[12]+pdata[13]+pdata[14])+pdata[15]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[8]-2.0f*(pdata[10]+pdata[11])+pdata[12]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[11])+pdata[13]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[10])+pdata[14]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[8]+pdata[11]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[10]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[9]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[8]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax); \
   value += (-10.0f*pdata[16]+4.0f*(pdata[17]+pdata[18]+pdata[19])-2.0f*(pdata[20]+pdata[21]+pdata[22])+pdata[23]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[16]-2.0f*(pdata[18]+pdata[19])+pdata[20]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[19])+pdata[21]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[18])+pdata[22]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[16]+pdata[19]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[18]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[17]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[16]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f); \
   value += (-10.0f*pdata[24]+4.0f*(pdata[25]+pdata[26]+pdata[27])-2.0f*(pdata[28]+pdata[29]+pdata[30])+pdata[31]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[24]-2.0f*(pdata[26]+pdata[27])+pdata[28]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[27])+pdata[29]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[26])+pdata[30]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[24]+pdata[27]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[26]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[25]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[24]) * rho22(alpha-1, beta-1, gamma-1);

#define CONV_PPIPED_OLD_INPLACE_RegionOne_Odd(value) \
   pdata[0] = Fetch(startX, startY, startZ); \
   pdata[1] = Fetch(startX+1, startY, startZ-1); \
   pdata[2] = Fetch(startX, startY+1, startZ-1); \
   pdata[3] = Fetch(startX, startY, startZ+1); \
   pdata[4] = Fetch(startX-1, startY, startZ); \
   pdata[5] = Fetch(startX, startY-1, startZ); \
   pdata[6] = Fetch(startX, startY, startZ-2); \
   pdata[7] = Fetch(startX, startY, startZ-1); \
   pdata[8] = Fetch(startX+1, startY+1, startZ+1); \
   pdata[9] = Fetch(startX, startY+1, startZ); \
   pdata[10] = Fetch(startX, startY, startZ+2); \
   pdata[11] = Fetch(startX+1, startY+1, startZ+2); \
   pdata[12] = Fetch(startX+1, startY+1, startZ+3); \
   pdata[13] = Fetch(startX+1, startY+2, startZ+1); \
   pdata[14] = Fetch(startX, startY+1, startZ+1); \
   pdata[15] = Fetch(startX, startY+1, startZ+2); \
   pdata[16] = Fetch(startX+1, startY+1, startZ-1); \
   pdata[17] = Fetch(startX+1, startY+1, startZ); \
   pdata[18] = Fetch(startX+1, startY, startZ-2); \
   pdata[19] = Fetch(startX, startY+1, startZ-2); \
   pdata[20] = Fetch(startX+1, startY+1, startZ-3); \
   pdata[21] = Fetch(startX+1, startY+2, startZ-1); \
   pdata[22] = Fetch(startX+2, startY+1, startZ-1); \
   pdata[23] = Fetch(startX+1, startY+1, startZ-2); \
   pdata[24] = Fetch(startX+1, startY, startZ); \
   pdata[25] = Fetch(startX+1, startY, startZ+1); \
   pdata[26] = Fetch(startX+2, startY+1, startZ+1); \
   pdata[27] = Fetch(startX+2, startY, startZ-1); \
   pdata[28] = Fetch(startX+2, startY, startZ); \
   pdata[29] = Fetch(startX+1, startY-1, startZ); \
   pdata[30] = Fetch(startX+1, startY, startZ+2); \
   pdata[31] = Fetch(startX+2, startY, startZ+1); \
   alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f; \
   value += (-10.0f*pdata[0]+4.0f*(pdata[1]+pdata[2]+pdata[3])-2.0f*(pdata[4]+pdata[5]+pdata[6])+pdata[7]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[0]-2.0f*(pdata[2]+pdata[3])+pdata[4]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[3])+pdata[5]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[2])+pdata[6]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[0]+pdata[3]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[2]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[1]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[0]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f; \
   value += (-10.0f*pdata[8]+4.0f*(pdata[9]+pdata[10]+pdata[11])-2.0f*(pdata[12]+pdata[13]+pdata[14])+pdata[15]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[8]-2.0f*(pdata[10]+pdata[11])+pdata[12]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[11])+pdata[13]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[10])+pdata[14]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[8]+pdata[11]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[10]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[9]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[8]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax); \
   value += (-10.0f*pdata[16]+4.0f*(pdata[17]+pdata[18]+pdata[19])-2.0f*(pdata[20]+pdata[21]+pdata[22])+pdata[23]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[16]-2.0f*(pdata[18]+pdata[19])+pdata[20]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[19])+pdata[21]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[18])+pdata[22]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[16]+pdata[19]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[18]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[17]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[16]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f); \
   value += (-10.0f*pdata[24]+4.0f*(pdata[25]+pdata[26]+pdata[27])-2.0f*(pdata[28]+pdata[29]+pdata[30])+pdata[31]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[24]-2.0f*(pdata[26]+pdata[27])+pdata[28]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[27])+pdata[29]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[26])+pdata[30]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[24]+pdata[27]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[26]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[25]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[24]) * rho22(alpha-1, beta-1, gamma-1);

#define CONV_PPIPED_OLD_INPLACE_RegionTwo_Even(value) \
   pdata[0] = Fetch(startX, startY, startZ); \
   pdata[1] = Fetch(startX-1, startY, startZ-1); \
   pdata[2] = Fetch(startX, startY-1, startZ-1); \
   pdata[3] = Fetch(startX-1, startY-1, startZ+1); \
   pdata[4] = Fetch(startX, startY-1, startZ); \
   pdata[5] = Fetch(startX-1, startY, startZ); \
   pdata[6] = Fetch(startX, startY, startZ-2); \
   pdata[7] = Fetch(startX-1, startY-1, startZ-1); \
   pdata[8] = Fetch(startX, startY, startZ+1); \
   pdata[9] = Fetch(startX+1, startY, startZ); \
   pdata[10] = Fetch(startX, startY, startZ+2); \
   pdata[11] = Fetch(startX+1, startY+1, startZ+2); \
   pdata[12] = Fetch(startX, startY, startZ+3); \
   pdata[13] = Fetch(startX+1, startY, startZ+1); \
   pdata[14] = Fetch(startX, startY-1, startZ+1); \
   pdata[15] = Fetch(startX+1, startY, startZ+2); \
   pdata[16] = Fetch(startX, startY, startZ-1); \
   pdata[17] = Fetch(startX+1, startY+1, startZ); \
   pdata[18] = Fetch(startX, startY+1, startZ-2); \
   pdata[19] = Fetch(startX+1, startY, startZ-2); \
   pdata[20] = Fetch(startX, startY, startZ-3); \
   pdata[21] = Fetch(startX+1, startY, startZ-1); \
   pdata[22] = Fetch(startX, startY+1, startZ-1); \
   pdata[23] = Fetch(startX+1, startY+1, startZ-2); \
   pdata[24] = Fetch(startX, startY+1, startZ); \
   pdata[25] = Fetch(startX-1, startY, startZ+1); \
   pdata[26] = Fetch(startX, startY+1, startZ+1); \
   pdata[27] = Fetch(startX-1, startY+1, startZ-1); \
   pdata[28] = Fetch(startX, startY+2, startZ); \
   pdata[29] = Fetch(startX-1, startY+1, startZ); \
   pdata[30] = Fetch(startX, startY+1, startZ+2); \
   pdata[31] = Fetch(startX-1, startY+1, startZ+1); \
   alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f; \
   value += (-10.0f*pdata[0]+4.0f*(pdata[1]+pdata[2]+pdata[3])-2.0f*(pdata[4]+pdata[5]+pdata[6])+pdata[7]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[0]-2.0f*(pdata[2]+pdata[3])+pdata[4]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[3])+pdata[5]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[2])+pdata[6]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[0]+pdata[3]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[2]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[1]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[0]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f; \
   value += (-10.0f*pdata[8]+4.0f*(pdata[9]+pdata[10]+pdata[11])-2.0f*(pdata[12]+pdata[13]+pdata[14])+pdata[15]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[8]-2.0f*(pdata[10]+pdata[11])+pdata[12]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[11])+pdata[13]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[10])+pdata[14]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[8]+pdata[11]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[10]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[9]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[8]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax); \
   value += (-10.0f*pdata[16]+4.0f*(pdata[17]+pdata[18]+pdata[19])-2.0f*(pdata[20]+pdata[21]+pdata[22])+pdata[23]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[16]-2.0f*(pdata[18]+pdata[19])+pdata[20]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[19])+pdata[21]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[18])+pdata[22]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[16]+pdata[19]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[18]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[17]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[16]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f); \
   value += (-10.0f*pdata[24]+4.0f*(pdata[25]+pdata[26]+pdata[27])-2.0f*(pdata[28]+pdata[29]+pdata[30])+pdata[31]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[24]-2.0f*(pdata[26]+pdata[27])+pdata[28]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[27])+pdata[29]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[26])+pdata[30]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[24]+pdata[27]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[26]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[25]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[24]) * rho22(alpha-1, beta-1, gamma-1);

#define CONV_PPIPED_OLD_INPLACE_RegionTwo_Odd(value) \
   pdata[0] = Fetch(startX, startY, startZ); \
   pdata[1] = Fetch(startX, startY+1, startZ-1); \
   pdata[2] = Fetch(startX+1, startY, startZ-1); \
   pdata[3] = Fetch(startX, startY, startZ+1); \
   pdata[4] = Fetch(startX, startY-1, startZ); \
   pdata[5] = Fetch(startX-1, startY, startZ); \
   pdata[6] = Fetch(startX, startY, startZ-2); \
   pdata[7] = Fetch(startX, startY, startZ-1); \
   pdata[8] = Fetch(startX+1, startY+1, startZ+1); \
   pdata[9] = Fetch(startX+1, startY, startZ); \
   pdata[10] = Fetch(startX, startY, startZ+2); \
   pdata[11] = Fetch(startX+1, startY+1, startZ+2); \
   pdata[12] = Fetch(startX+1, startY+1, startZ+3); \
   pdata[13] = Fetch(startX+2, startY+1, startZ+1); \
   pdata[14] = Fetch(startX+1, startY, startZ+1); \
   pdata[15] = Fetch(startX+1, startY, startZ+2); \
   pdata[16] = Fetch(startX+1, startY+1, startZ-1); \
   pdata[17] = Fetch(startX+1, startY+1, startZ); \
   pdata[18] = Fetch(startX, startY+1, startZ-2); \
   pdata[19] = Fetch(startX+1, startY, startZ-2); \
   pdata[20] = Fetch(startX+1, startY+1, startZ-3); \
   pdata[21] = Fetch(startX+2, startY+1, startZ-1); \
   pdata[22] = Fetch(startX+1, startY+2, startZ-1); \
   pdata[23] = Fetch(startX+1, startY+1, startZ-2); \
   pdata[24] = Fetch(startX, startY+1, startZ); \
   pdata[25] = Fetch(startX, startY+1, startZ+1); \
   pdata[26] = Fetch(startX+1, startY+2, startZ+1); \
   pdata[27] = Fetch(startX, startY+2, startZ-1); \
   pdata[28] = Fetch(startX, startY+2, startZ); \
   pdata[29] = Fetch(startX-1, startY+1, startZ); \
   pdata[30] = Fetch(startX, startY+1, startZ+2); \
   pdata[31] = Fetch(startX, startY+2, startZ+1); \
   alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f; \
   value += (-10.0f*pdata[0]+4.0f*(pdata[1]+pdata[2]+pdata[3])-2.0f*(pdata[4]+pdata[5]+pdata[6])+pdata[7]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[0]-2.0f*(pdata[2]+pdata[3])+pdata[4]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[3])+pdata[5]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[2])+pdata[6]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[0]+pdata[3]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[2]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[1]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[0]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f; \
   value += (-10.0f*pdata[8]+4.0f*(pdata[9]+pdata[10]+pdata[11])-2.0f*(pdata[12]+pdata[13]+pdata[14])+pdata[15]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[8]-2.0f*(pdata[10]+pdata[11])+pdata[12]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[11])+pdata[13]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[10])+pdata[14]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[8]+pdata[11]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[10]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[9]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[8]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax); \
   value += (-10.0f*pdata[16]+4.0f*(pdata[17]+pdata[18]+pdata[19])-2.0f*(pdata[20]+pdata[21]+pdata[22])+pdata[23]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[16]-2.0f*(pdata[18]+pdata[19])+pdata[20]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[19])+pdata[21]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[18])+pdata[22]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[16]+pdata[19]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[18]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[17]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[16]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f); \
   value += (-10.0f*pdata[24]+4.0f*(pdata[25]+pdata[26]+pdata[27])-2.0f*(pdata[28]+pdata[29]+pdata[30])+pdata[31]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[24]-2.0f*(pdata[26]+pdata[27])+pdata[28]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[27])+pdata[29]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[26])+pdata[30]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[24]+pdata[27]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[26]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[25]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[24]) * rho22(alpha-1, beta-1, gamma-1);

#define CONV_PPIPED_OLD_INPLACE_RegionThree_Even(value) \
   pdata[0] = Fetch(startX, startY, startZ); \
   pdata[1] = Fetch(startX-1, startY, startZ-1); \
   pdata[2] = Fetch(startX-1, startY-1, startZ+1); \
   pdata[3] = Fetch(startX, startY-1, startZ-1); \
   pdata[4] = Fetch(startX, startY-1, startZ); \
   pdata[5] = Fetch(startX, startY, startZ-2); \
   pdata[6] = Fetch(startX-1, startY, startZ); \
   pdata[7] = Fetch(startX-1, startY-1, startZ-1); \
   pdata[8] = Fetch(startX, startY, startZ+1); \
   pdata[9] = Fetch(startX, startY, startZ+2); \
   pdata[10] = Fetch(startX+1, startY, startZ); \
   pdata[11] = Fetch(startX+1, startY+1, startZ+2); \
   pdata[12] = Fetch(startX+1, startY, startZ+1); \
   pdata[13] = Fetch(startX, startY, startZ+3); \
   pdata[14] = Fetch(startX, startY-1, startZ+1); \
   pdata[15] = Fetch(startX+1, startY, startZ+2); \
   pdata[16] = Fetch(startX-1, startY, startZ+1); \
   pdata[17] = Fetch(startX, startY+1, startZ+2); \
   pdata[18] = Fetch(startX-1, startY+1, startZ); \
   pdata[19] = Fetch(startX-1, startY, startZ+2); \
   pdata[20] = Fetch(startX-2, startY, startZ+1); \
   pdata[21] = Fetch(startX-1, startY, startZ+3); \
   pdata[22] = Fetch(startX-1, startY+1, startZ+1); \
   pdata[23] = Fetch(startX-1, startY+1, startZ+2); \
   pdata[24] = Fetch(startX, startY+1, startZ); \
   pdata[25] = Fetch(startX, startY, startZ-1); \
   pdata[26] = Fetch(startX, startY+1, startZ+1); \
   pdata[27] = Fetch(startX-1, startY+1, startZ-1); \
   pdata[28] = Fetch(startX, startY+2, startZ); \
   pdata[29] = Fetch(startX, startY+1, startZ-2); \
   pdata[30] = Fetch(startX+1, startY+1, startZ); \
   pdata[31] = Fetch(startX, startY+1, startZ-1); \
   alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f; \
   value += (-10.0f*pdata[0]+4.0f*(pdata[1]+pdata[2]+pdata[3])-2.0f*(pdata[4]+pdata[5]+pdata[6])+pdata[7]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[0]-2.0f*(pdata[2]+pdata[3])+pdata[4]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[3])+pdata[5]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[2])+pdata[6]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[0]+pdata[3]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[2]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[1]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[0]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f; \
   value += (-10.0f*pdata[8]+4.0f*(pdata[9]+pdata[10]+pdata[11])-2.0f*(pdata[12]+pdata[13]+pdata[14])+pdata[15]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[8]-2.0f*(pdata[10]+pdata[11])+pdata[12]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[11])+pdata[13]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[10])+pdata[14]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[8]+pdata[11]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[10]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[9]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[8]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax); \
   value += (-10.0f*pdata[16]+4.0f*(pdata[17]+pdata[18]+pdata[19])-2.0f*(pdata[20]+pdata[21]+pdata[22])+pdata[23]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[16]-2.0f*(pdata[18]+pdata[19])+pdata[20]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[19])+pdata[21]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[18])+pdata[22]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[16]+pdata[19]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[18]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[17]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[16]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f); \
   value += (-10.0f*pdata[24]+4.0f*(pdata[25]+pdata[26]+pdata[27])-2.0f*(pdata[28]+pdata[29]+pdata[30])+pdata[31]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[24]-2.0f*(pdata[26]+pdata[27])+pdata[28]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[27])+pdata[29]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[26])+pdata[30]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[24]+pdata[27]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[26]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[25]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[24]) * rho22(alpha-1, beta-1, gamma-1);

#define CONV_PPIPED_OLD_INPLACE_RegionThree_Odd(value) \
   pdata[0] = Fetch(startX, startY, startZ); \
   pdata[1] = Fetch(startX, startY+1, startZ-1); \
   pdata[2] = Fetch(startX, startY, startZ+1); \
   pdata[3] = Fetch(startX+1, startY, startZ-1); \
   pdata[4] = Fetch(startX, startY-1, startZ); \
   pdata[5] = Fetch(startX, startY, startZ-2); \
   pdata[6] = Fetch(startX-1, startY, startZ); \
   pdata[7] = Fetch(startX, startY, startZ-1); \
   pdata[8] = Fetch(startX+1, startY+1, startZ+1); \
   pdata[9] = Fetch(startX, startY, startZ+2); \
   pdata[10] = Fetch(startX+1, startY, startZ); \
   pdata[11] = Fetch(startX+1, startY+1, startZ+2); \
   pdata[12] = Fetch(startX+2, startY+1, startZ+1); \
   pdata[13] = Fetch(startX+1, startY+1, startZ+3); \
   pdata[14] = Fetch(startX+1, startY, startZ+1); \
   pdata[15] = Fetch(startX+1, startY, startZ+2); \
   pdata[16] = Fetch(startX, startY+1, startZ+1); \
   pdata[17] = Fetch(startX, startY+1, startZ+2); \
   pdata[18] = Fetch(startX-1, startY+1, startZ); \
   pdata[19] = Fetch(startX-1, startY, startZ+2); \
   pdata[20] = Fetch(startX-1, startY+1, startZ+1); \
   pdata[21] = Fetch(startX, startY+1, startZ+3); \
   pdata[22] = Fetch(startX, startY+2, startZ+1); \
   pdata[23] = Fetch(startX-1, startY+1, startZ+2); \
   pdata[24] = Fetch(startX, startY+1, startZ); \
   pdata[25] = Fetch(startX+1, startY+1, startZ-1); \
   pdata[26] = Fetch(startX+1, startY+2, startZ+1); \
   pdata[27] = Fetch(startX, startY+2, startZ-1); \
   pdata[28] = Fetch(startX, startY+2, startZ); \
   pdata[29] = Fetch(startX, startY+1, startZ-2); \
   pdata[30] = Fetch(startX+1, startY+1, startZ); \
   pdata[31] = Fetch(startX+1, startY+2, startZ-1); \
   alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f; \
   value += (-10.0f*pdata[0]+4.0f*(pdata[1]+pdata[2]+pdata[3])-2.0f*(pdata[4]+pdata[5]+pdata[6])+pdata[7]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[0]-2.0f*(pdata[2]+pdata[3])+pdata[4]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[3])+pdata[5]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[2])+pdata[6]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[0]+pdata[3]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[2]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[1]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[0]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f; \
   value += (-10.0f*pdata[8]+4.0f*(pdata[9]+pdata[10]+pdata[11])-2.0f*(pdata[12]+pdata[13]+pdata[14])+pdata[15]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[8]-2.0f*(pdata[10]+pdata[11])+pdata[12]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[11])+pdata[13]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[10])+pdata[14]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[8]+pdata[11]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[10]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[9]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[8]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax); \
   value += (-10.0f*pdata[16]+4.0f*(pdata[17]+pdata[18]+pdata[19])-2.0f*(pdata[20]+pdata[21]+pdata[22])+pdata[23]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[16]-2.0f*(pdata[18]+pdata[19])+pdata[20]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[19])+pdata[21]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[18])+pdata[22]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[16]+pdata[19]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[18]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[17]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[16]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f); \
   value += (-10.0f*pdata[24]+4.0f*(pdata[25]+pdata[26]+pdata[27])-2.0f*(pdata[28]+pdata[29]+pdata[30])+pdata[31]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[24]-2.0f*(pdata[26]+pdata[27])+pdata[28]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[27])+pdata[29]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[26])+pdata[30]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[24]+pdata[27]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[26]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[25]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[24]) * rho22(alpha-1, beta-1, gamma-1);

#define CONV_PPIPED_OLD_INPLACE_RegionFour_Even(value) \
   pdata[0] = Fetch(startX, startY, startZ); \
   pdata[1] = Fetch(startX, startY-1, startZ-1); \
   pdata[2] = Fetch(startX-1, startY-1, startZ+1); \
   pdata[3] = Fetch(startX-1, startY, startZ-1); \
   pdata[4] = Fetch(startX-1, startY, startZ); \
   pdata[5] = Fetch(startX, startY, startZ-2); \
   pdata[6] = Fetch(startX, startY-1, startZ); \
   pdata[7] = Fetch(startX-1, startY-1, startZ-1); \
   pdata[8] = Fetch(startX, startY, startZ+1); \
   pdata[9] = Fetch(startX, startY, startZ+2); \
   pdata[10] = Fetch(startX, startY+1, startZ); \
   pdata[11] = Fetch(startX+1, startY+1, startZ+2); \
   pdata[12] = Fetch(startX, startY+1, startZ+1); \
   pdata[13] = Fetch(startX, startY, startZ+3); \
   pdata[14] = Fetch(startX-1, startY, startZ+1); \
   pdata[15] = Fetch(startX, startY+1, startZ+2); \
   pdata[16] = Fetch(startX, startY-1, startZ+1); \
   pdata[17] = Fetch(startX+1, startY, startZ+2); \
   pdata[18] = Fetch(startX+1, startY-1, startZ); \
   pdata[19] = Fetch(startX, startY-1, startZ+2); \
   pdata[20] = Fetch(startX, startY-2, startZ+1); \
   pdata[21] = Fetch(startX, startY-1, startZ+3); \
   pdata[22] = Fetch(startX+1, startY-1, startZ+1); \
   pdata[23] = Fetch(startX+1, startY-1, startZ+2); \
   pdata[24] = Fetch(startX+1, startY, startZ); \
   pdata[25] = Fetch(startX, startY, startZ-1); \
   pdata[26] = Fetch(startX+1, startY, startZ+1); \
   pdata[27] = Fetch(startX+1, startY-1, startZ-1); \
   pdata[28] = Fetch(startX+2, startY, startZ); \
   pdata[29] = Fetch(startX+1, startY, startZ-2); \
   pdata[30] = Fetch(startX+1, startY+1, startZ); \
   pdata[31] = Fetch(startX+1, startY, startZ-1); \
   alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f; \
   value += (-10.0f*pdata[0]+4.0f*(pdata[1]+pdata[2]+pdata[3])-2.0f*(pdata[4]+pdata[5]+pdata[6])+pdata[7]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[0]-2.0f*(pdata[2]+pdata[3])+pdata[4]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[3])+pdata[5]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[2])+pdata[6]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[0]+pdata[3]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[2]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[1]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[0]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f; \
   value += (-10.0f*pdata[8]+4.0f*(pdata[9]+pdata[10]+pdata[11])-2.0f*(pdata[12]+pdata[13]+pdata[14])+pdata[15]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[8]-2.0f*(pdata[10]+pdata[11])+pdata[12]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[11])+pdata[13]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[10])+pdata[14]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[8]+pdata[11]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[10]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[9]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[8]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax); \
   value += (-10.0f*pdata[16]+4.0f*(pdata[17]+pdata[18]+pdata[19])-2.0f*(pdata[20]+pdata[21]+pdata[22])+pdata[23]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[16]-2.0f*(pdata[18]+pdata[19])+pdata[20]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[19])+pdata[21]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[18])+pdata[22]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[16]+pdata[19]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[18]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[17]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[16]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f); \
   value += (-10.0f*pdata[24]+4.0f*(pdata[25]+pdata[26]+pdata[27])-2.0f*(pdata[28]+pdata[29]+pdata[30])+pdata[31]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[24]-2.0f*(pdata[26]+pdata[27])+pdata[28]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[27])+pdata[29]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[26])+pdata[30]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[24]+pdata[27]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[26]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[25]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[24]) * rho22(alpha-1, beta-1, gamma-1);

#define CONV_PPIPED_OLD_INPLACE_RegionFour_Odd(value) \
   pdata[0] = Fetch(startX, startY, startZ); \
   pdata[1] = Fetch(startX+1, startY, startZ-1); \
   pdata[2] = Fetch(startX, startY, startZ+1); \
   pdata[3] = Fetch(startX, startY+1, startZ-1); \
   pdata[4] = Fetch(startX-1, startY, startZ); \
   pdata[5] = Fetch(startX, startY, startZ-2); \
   pdata[6] = Fetch(startX, startY-1, startZ); \
   pdata[7] = Fetch(startX, startY, startZ-1); \
   pdata[8] = Fetch(startX+1, startY+1, startZ+1); \
   pdata[9] = Fetch(startX, startY, startZ+2); \
   pdata[10] = Fetch(startX, startY+1, startZ); \
   pdata[11] = Fetch(startX+1, startY+1, startZ+2); \
   pdata[12] = Fetch(startX+1, startY+2, startZ+1); \
   pdata[13] = Fetch(startX+1, startY+1, startZ+3); \
   pdata[14] = Fetch(startX, startY+1, startZ+1); \
   pdata[15] = Fetch(startX, startY+1, startZ+2); \
   pdata[16] = Fetch(startX+1, startY, startZ+1); \
   pdata[17] = Fetch(startX+1, startY, startZ+2); \
   pdata[18] = Fetch(startX+1, startY-1, startZ); \
   pdata[19] = Fetch(startX, startY-1, startZ+2); \
   pdata[20] = Fetch(startX+1, startY-1, startZ+1); \
   pdata[21] = Fetch(startX+1, startY, startZ+3); \
   pdata[22] = Fetch(startX+2, startY, startZ+1); \
   pdata[23] = Fetch(startX+1, startY-1, startZ+2); \
   pdata[24] = Fetch(startX+1, startY, startZ); \
   pdata[25] = Fetch(startX+1, startY+1, startZ-1); \
   pdata[26] = Fetch(startX+2, startY+1, startZ+1); \
   pdata[27] = Fetch(startX+2, startY, startZ-1); \
   pdata[28] = Fetch(startX+2, startY, startZ); \
   pdata[29] = Fetch(startX+1, startY, startZ-2); \
   pdata[30] = Fetch(startX+1, startY+1, startZ); \
   pdata[31] = Fetch(startX+2, startY+1, startZ-1); \
   alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f; \
   value += (-10.0f*pdata[0]+4.0f*(pdata[1]+pdata[2]+pdata[3])-2.0f*(pdata[4]+pdata[5]+pdata[6])+pdata[7]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[0]-2.0f*(pdata[2]+pdata[3])+pdata[4]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[3])+pdata[5]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[2])+pdata[6]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[0]+pdata[3]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[2]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[1]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[0]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f; \
   value += (-10.0f*pdata[8]+4.0f*(pdata[9]+pdata[10]+pdata[11])-2.0f*(pdata[12]+pdata[13]+pdata[14])+pdata[15]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[8]-2.0f*(pdata[10]+pdata[11])+pdata[12]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[11])+pdata[13]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[10])+pdata[14]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[8]+pdata[11]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[10]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[9]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[8]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax); \
   value += (-10.0f*pdata[16]+4.0f*(pdata[17]+pdata[18]+pdata[19])-2.0f*(pdata[20]+pdata[21]+pdata[22])+pdata[23]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[16]-2.0f*(pdata[18]+pdata[19])+pdata[20]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[19])+pdata[21]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[18])+pdata[22]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[16]+pdata[19]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[18]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[17]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[16]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f); \
   value += (-10.0f*pdata[24]+4.0f*(pdata[25]+pdata[26]+pdata[27])-2.0f*(pdata[28]+pdata[29]+pdata[30])+pdata[31]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[24]-2.0f*(pdata[26]+pdata[27])+pdata[28]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[27])+pdata[29]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[26])+pdata[30]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[24]+pdata[27]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[26]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[25]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[24]) * rho22(alpha-1, beta-1, gamma-1);

#define CONV_PPIPED_OLD_INPLACE_RegionFive_Even(value) \
   pdata[0] = Fetch(startX, startY, startZ); \
   pdata[1] = Fetch(startX-1, startY-1, startZ+1); \
   pdata[2] = Fetch(startX, startY-1, startZ-1); \
   pdata[3] = Fetch(startX-1, startY, startZ-1); \
   pdata[4] = Fetch(startX, startY, startZ-2); \
   pdata[5] = Fetch(startX-1, startY, startZ); \
   pdata[6] = Fetch(startX, startY-1, startZ); \
   pdata[7] = Fetch(startX-1, startY-1, startZ-1); \
   pdata[8] = Fetch(startX, startY, startZ+1); \
   pdata[9] = Fetch(startX+1, startY, startZ); \
   pdata[10] = Fetch(startX, startY+1, startZ); \
   pdata[11] = Fetch(startX+1, startY+1, startZ+2); \
   pdata[12] = Fetch(startX, startY+1, startZ+1); \
   pdata[13] = Fetch(startX+1, startY, startZ+1); \
   pdata[14] = Fetch(startX, startY, startZ-1); \
   pdata[15] = Fetch(startX+1, startY+1, startZ); \
   pdata[16] = Fetch(startX, startY-1, startZ+1); \
   pdata[17] = Fetch(startX+1, startY, startZ+2); \
   pdata[18] = Fetch(startX, startY-1, startZ+2); \
   pdata[19] = Fetch(startX+1, startY-1, startZ); \
   pdata[20] = Fetch(startX, startY-2, startZ+1); \
   pdata[21] = Fetch(startX+1, startY-1, startZ+1); \
   pdata[22] = Fetch(startX, startY-1, startZ+3); \
   pdata[23] = Fetch(startX+1, startY-1, startZ+2); \
   pdata[24] = Fetch(startX, startY, startZ+2); \
   pdata[25] = Fetch(startX-1, startY, startZ+1); \
   pdata[26] = Fetch(startX, startY, startZ+3); \
   pdata[27] = Fetch(startX-1, startY-1, startZ+3); \
   pdata[28] = Fetch(startX, startY, startZ+4); \
   pdata[29] = Fetch(startX-1, startY, startZ+2); \
   pdata[30] = Fetch(startX, startY+1, startZ+2); \
   pdata[31] = Fetch(startX-1, startY, startZ+3); \
   alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f; \
   value += (-10.0f*pdata[0]+4.0f*(pdata[1]+pdata[2]+pdata[3])-2.0f*(pdata[4]+pdata[5]+pdata[6])+pdata[7]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[0]-2.0f*(pdata[2]+pdata[3])+pdata[4]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[3])+pdata[5]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[2])+pdata[6]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[0]+pdata[3]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[2]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[1]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[0]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f; \
   value += (-10.0f*pdata[8]+4.0f*(pdata[9]+pdata[10]+pdata[11])-2.0f*(pdata[12]+pdata[13]+pdata[14])+pdata[15]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[8]-2.0f*(pdata[10]+pdata[11])+pdata[12]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[11])+pdata[13]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[10])+pdata[14]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[8]+pdata[11]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[10]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[9]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[8]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax); \
   value += (-10.0f*pdata[16]+4.0f*(pdata[17]+pdata[18]+pdata[19])-2.0f*(pdata[20]+pdata[21]+pdata[22])+pdata[23]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[16]-2.0f*(pdata[18]+pdata[19])+pdata[20]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[19])+pdata[21]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[18])+pdata[22]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[16]+pdata[19]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[18]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[17]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[16]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f); \
   value += (-10.0f*pdata[24]+4.0f*(pdata[25]+pdata[26]+pdata[27])-2.0f*(pdata[28]+pdata[29]+pdata[30])+pdata[31]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[24]-2.0f*(pdata[26]+pdata[27])+pdata[28]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[27])+pdata[29]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[26])+pdata[30]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[24]+pdata[27]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[26]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[25]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[24]) * rho22(alpha-1, beta-1, gamma-1);

#define CONV_PPIPED_OLD_INPLACE_RegionFive_Odd(value) \
   pdata[0] = Fetch(startX, startY, startZ); \
   pdata[1] = Fetch(startX, startY, startZ+1); \
   pdata[2] = Fetch(startX+1, startY, startZ-1); \
   pdata[3] = Fetch(startX, startY+1, startZ-1); \
   pdata[4] = Fetch(startX, startY, startZ-2); \
   pdata[5] = Fetch(startX-1, startY, startZ); \
   pdata[6] = Fetch(startX, startY-1, startZ); \
   pdata[7] = Fetch(startX, startY, startZ-1); \
   pdata[8] = Fetch(startX+1, startY+1, startZ+1); \
   pdata[9] = Fetch(startX+1, startY, startZ); \
   pdata[10] = Fetch(startX, startY+1, startZ); \
   pdata[11] = Fetch(startX+1, startY+1, startZ+2); \
   pdata[12] = Fetch(startX+1, startY+2, startZ+1); \
   pdata[13] = Fetch(startX+2, startY+1, startZ+1); \
   pdata[14] = Fetch(startX+1, startY+1, startZ-1); \
   pdata[15] = Fetch(startX+1, startY+1, startZ); \
   pdata[16] = Fetch(startX+1, startY, startZ+1); \
   pdata[17] = Fetch(startX+1, startY, startZ+2); \
   pdata[18] = Fetch(startX, startY-1, startZ+2); \
   pdata[19] = Fetch(startX+1, startY-1, startZ); \
   pdata[20] = Fetch(startX+1, startY-1, startZ+1); \
   pdata[21] = Fetch(startX+2, startY, startZ+1); \
   pdata[22] = Fetch(startX+1, startY, startZ+3); \
   pdata[23] = Fetch(startX+1, startY-1, startZ+2); \
   pdata[24] = Fetch(startX, startY, startZ+2); \
   pdata[25] = Fetch(startX, startY+1, startZ+1); \
   pdata[26] = Fetch(startX+1, startY+1, startZ+3); \
   pdata[27] = Fetch(startX, startY, startZ+3); \
   pdata[28] = Fetch(startX, startY, startZ+4); \
   pdata[29] = Fetch(startX-1, startY, startZ+2); \
   pdata[30] = Fetch(startX, startY+1, startZ+2); \
   pdata[31] = Fetch(startX, startY+1, startZ+3); \
   alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f; \
   value += (-10.0f*pdata[0]+4.0f*(pdata[1]+pdata[2]+pdata[3])-2.0f*(pdata[4]+pdata[5]+pdata[6])+pdata[7]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[0]-2.0f*(pdata[2]+pdata[3])+pdata[4]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[3])+pdata[5]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[2])+pdata[6]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[0]+pdata[3]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[2]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[1]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[0]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f; \
   value += (-10.0f*pdata[8]+4.0f*(pdata[9]+pdata[10]+pdata[11])-2.0f*(pdata[12]+pdata[13]+pdata[14])+pdata[15]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[8]-2.0f*(pdata[10]+pdata[11])+pdata[12]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[11])+pdata[13]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[10])+pdata[14]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[8]+pdata[11]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[10]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[9]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[8]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax); \
   value += (-10.0f*pdata[16]+4.0f*(pdata[17]+pdata[18]+pdata[19])-2.0f*(pdata[20]+pdata[21]+pdata[22])+pdata[23]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[16]-2.0f*(pdata[18]+pdata[19])+pdata[20]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[19])+pdata[21]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[18])+pdata[22]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[16]+pdata[19]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[18]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[17]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[16]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f); \
   value += (-10.0f*pdata[24]+4.0f*(pdata[25]+pdata[26]+pdata[27])-2.0f*(pdata[28]+pdata[29]+pdata[30])+pdata[31]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[24]-2.0f*(pdata[26]+pdata[27])+pdata[28]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[27])+pdata[29]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[26])+pdata[30]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[24]+pdata[27]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[26]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[25]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[24]) * rho22(alpha-1, beta-1, gamma-1);

#define CONV_PPIPED_OLD_INPLACE_RegionSix_Even(value) \
   pdata[0] = Fetch(startX, startY, startZ); \
   pdata[1] = Fetch(startX-1, startY-1, startZ+1); \
   pdata[2] = Fetch(startX-1, startY, startZ-1); \
   pdata[3] = Fetch(startX, startY-1, startZ-1); \
   pdata[4] = Fetch(startX, startY, startZ-2); \
   pdata[5] = Fetch(startX, startY-1, startZ); \
   pdata[6] = Fetch(startX-1, startY, startZ); \
   pdata[7] = Fetch(startX-1, startY-1, startZ-1); \
   pdata[8] = Fetch(startX, startY, startZ+1); \
   pdata[9] = Fetch(startX, startY+1, startZ); \
   pdata[10] = Fetch(startX+1, startY, startZ); \
   pdata[11] = Fetch(startX+1, startY+1, startZ+2); \
   pdata[12] = Fetch(startX+1, startY, startZ+1); \
   pdata[13] = Fetch(startX, startY+1, startZ+1); \
   pdata[14] = Fetch(startX, startY, startZ-1); \
   pdata[15] = Fetch(startX+1, startY+1, startZ); \
   pdata[16] = Fetch(startX-1, startY, startZ+1); \
   pdata[17] = Fetch(startX, startY+1, startZ+2); \
   pdata[18] = Fetch(startX-1, startY, startZ+2); \
   pdata[19] = Fetch(startX-1, startY+1, startZ); \
   pdata[20] = Fetch(startX-2, startY, startZ+1); \
   pdata[21] = Fetch(startX-1, startY+1, startZ+1); \
   pdata[22] = Fetch(startX-1, startY, startZ+3); \
   pdata[23] = Fetch(startX-1, startY+1, startZ+2); \
   pdata[24] = Fetch(startX, startY, startZ+2); \
   pdata[25] = Fetch(startX, startY-1, startZ+1); \
   pdata[26] = Fetch(startX, startY, startZ+3); \
   pdata[27] = Fetch(startX-1, startY-1, startZ+3); \
   pdata[28] = Fetch(startX, startY, startZ+4); \
   pdata[29] = Fetch(startX, startY-1, startZ+2); \
   pdata[30] = Fetch(startX+1, startY, startZ+2); \
   pdata[31] = Fetch(startX, startY-1, startZ+3); \
   alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f; \
   value += (-10.0f*pdata[0]+4.0f*(pdata[1]+pdata[2]+pdata[3])-2.0f*(pdata[4]+pdata[5]+pdata[6])+pdata[7]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[0]-2.0f*(pdata[2]+pdata[3])+pdata[4]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[3])+pdata[5]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[2])+pdata[6]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[0]+pdata[3]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[2]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[1]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[0]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f; \
   value += (-10.0f*pdata[8]+4.0f*(pdata[9]+pdata[10]+pdata[11])-2.0f*(pdata[12]+pdata[13]+pdata[14])+pdata[15]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[8]-2.0f*(pdata[10]+pdata[11])+pdata[12]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[11])+pdata[13]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[10])+pdata[14]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[8]+pdata[11]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[10]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[9]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[8]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax); \
   value += (-10.0f*pdata[16]+4.0f*(pdata[17]+pdata[18]+pdata[19])-2.0f*(pdata[20]+pdata[21]+pdata[22])+pdata[23]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[16]-2.0f*(pdata[18]+pdata[19])+pdata[20]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[19])+pdata[21]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[18])+pdata[22]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[16]+pdata[19]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[18]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[17]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[16]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f); \
   value += (-10.0f*pdata[24]+4.0f*(pdata[25]+pdata[26]+pdata[27])-2.0f*(pdata[28]+pdata[29]+pdata[30])+pdata[31]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[24]-2.0f*(pdata[26]+pdata[27])+pdata[28]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[27])+pdata[29]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[26])+pdata[30]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[24]+pdata[27]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[26]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[25]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[24]) * rho22(alpha-1, beta-1, gamma-1);

#define CONV_PPIPED_OLD_INPLACE_RegionSix_Odd(value) \
   pdata[0] = Fetch(startX, startY, startZ); \
   pdata[1] = Fetch(startX, startY, startZ+1); \
   pdata[2] = Fetch(startX, startY+1, startZ-1); \
   pdata[3] = Fetch(startX+1, startY, startZ-1); \
   pdata[4] = Fetch(startX, startY, startZ-2); \
   pdata[5] = Fetch(startX, startY-1, startZ); \
   pdata[6] = Fetch(startX-1, startY, startZ); \
   pdata[7] = Fetch(startX, startY, startZ-1); \
   pdata[8] = Fetch(startX+1, startY+1, startZ+1); \
   pdata[9] = Fetch(startX, startY+1, startZ); \
   pdata[10] = Fetch(startX+1, startY, startZ); \
   pdata[11] = Fetch(startX+1, startY+1, startZ+2); \
   pdata[12] = Fetch(startX+2, startY+1, startZ+1); \
   pdata[13] = Fetch(startX+1, startY+2, startZ+1); \
   pdata[14] = Fetch(startX+1, startY+1, startZ-1); \
   pdata[15] = Fetch(startX+1, startY+1, startZ); \
   pdata[16] = Fetch(startX, startY+1, startZ+1); \
   pdata[17] = Fetch(startX, startY+1, startZ+2); \
   pdata[18] = Fetch(startX-1, startY, startZ+2); \
   pdata[19] = Fetch(startX-1, startY+1, startZ); \
   pdata[20] = Fetch(startX-1, startY+1, startZ+1); \
   pdata[21] = Fetch(startX, startY+2, startZ+1); \
   pdata[22] = Fetch(startX, startY+1, startZ+3); \
   pdata[23] = Fetch(startX-1, startY+1, startZ+2); \
   pdata[24] = Fetch(startX, startY, startZ+2); \
   pdata[25] = Fetch(startX+1, startY, startZ+1); \
   pdata[26] = Fetch(startX+1, startY+1, startZ+3); \
   pdata[27] = Fetch(startX, startY, startZ+3); \
   pdata[28] = Fetch(startX, startY, startZ+4); \
   pdata[29] = Fetch(startX, startY-1, startZ+2); \
   pdata[30] = Fetch(startX+1, startY, startZ+2); \
   pdata[31] = Fetch(startX+1, startY, startZ+3); \
   alpha=mymax-1.0f;beta=mymid-1.0f;gamma=mymin-1.0f; \
   value += (-10.0f*pdata[0]+4.0f*(pdata[1]+pdata[2]+pdata[3])-2.0f*(pdata[4]+pdata[5]+pdata[6])+pdata[7]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[0]-2.0f*(pdata[2]+pdata[3])+pdata[4]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[3])+pdata[5]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[0]-2.0f*(pdata[1]+pdata[2])+pdata[6]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[0]+pdata[3]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[2]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[0]+pdata[1]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[0]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=-mymin;beta=mymax-mymin-1.0f;gamma=mymid-mymin-1.0f; \
   value += (-10.0f*pdata[8]+4.0f*(pdata[9]+pdata[10]+pdata[11])-2.0f*(pdata[12]+pdata[13]+pdata[14])+pdata[15]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[8]-2.0f*(pdata[10]+pdata[11])+pdata[12]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[11])+pdata[13]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[8]-2.0f*(pdata[9]+pdata[10])+pdata[14]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[8]+pdata[11]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[10]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[8]+pdata[9]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[8]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymax+mymid);beta=(-mymax+mymin);gamma=(-mymax); \
   value += (-10.0f*pdata[16]+4.0f*(pdata[17]+pdata[18]+pdata[19])-2.0f*(pdata[20]+pdata[21]+pdata[22])+pdata[23]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[16]-2.0f*(pdata[18]+pdata[19])+pdata[20]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[19])+pdata[21]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[16]-2.0f*(pdata[17]+pdata[18])+pdata[22]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[16]+pdata[19]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[18]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[16]+pdata[17]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[16]) * rho22(alpha-1, beta-1, gamma-1);\
   alpha=(-mymid+mymin);beta=(-mymid);gamma=(mymax-mymid-1.0f); \
   value += (-10.0f*pdata[24]+4.0f*(pdata[25]+pdata[26]+pdata[27])-2.0f*(pdata[28]+pdata[29]+pdata[30])+pdata[31]) \
   * rho22(alpha, beta, gamma) + \
   (4.0f*pdata[24]-2.0f*(pdata[26]+pdata[27])+pdata[28]) \
   * rho22(alpha, beta, gamma-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[27])+pdata[29]) \
   * rho22(alpha, gamma, beta-1) + \
   (4.0f*pdata[24]-2.0f*(pdata[25]+pdata[26])+pdata[30]) \
   * rho22(beta, gamma, alpha-1) + \
   (-2.0f*pdata[24]+pdata[27]) * rho22(alpha, beta-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[26]) * rho22(beta, alpha-1, gamma-1) + \
   (-2.0f*pdata[24]+pdata[25]) * rho22(gamma, alpha-1, beta-1) + \
   (pdata[24]) * rho22(alpha-1, beta-1, gamma-1);

