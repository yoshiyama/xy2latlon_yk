#include <iostream>
#include <cmath>
#include <vector>
#include <utility>

std::pair<double, double> calc_lat_lon(double x, double y, double phi0_deg, double lambda0_deg) {
    // 平面直角座標系原点をラジアンに直す
    double phi0_rad = phi0_deg * M_PI / 180.0;
    double lambda0_rad = lambda0_deg * M_PI / 180.0;

    // 補助関数
    auto A_array = [](double n) -> std::vector<double> {
        double A0 = 1 + (n * n) / 4.0 + (n * n * n * n) / 64.0;
        double A1 = -(3.0 / 2.0) * (n - (n * n * n) / 8.0 - (n * n * n * n * n) / 64.0);
        double A2 = (15.0 / 16.0) * (n * n - (n * n * n * n) / 4.0);
        double A3 = -(35.0 / 48.0) * (n * n * n - (5.0 / 16.0) * (n * n * n * n * n));
        double A4 = (315.0 / 512.0) * (n * n * n * n);
        double A5 = -(693.0 / 1280.0) * (n * n * n * n * n);
        return std::vector<double>{A0, A1, A2, A3, A4, A5};
    };

    auto beta_array = [](double n) -> std::vector<double> {
        double b0 = std::nan(""); // dummy
        double b1 = (1.0 / 2.0) * n - (2.0 / 3.0) * (n * n) + (37.0 / 96.0) * (n * n * n) - (1.0 / 360.0) * (n * n * n * n) - (81.0 / 512.0) * (n * n * n * n * n);
        double b2 = (1.0 / 48.0) * (n * n) + (1.0 / 15.0) * (n * n * n) - (437.0 / 1440.0) * (n * n * n * n) + (46.0 / 105.0) * (n * n * n * n * n);
        double b3 = (17.0 / 480.0) * (n * n * n) - (37.0 / 840.0) * (n * n * n * n) - (209.0 / 4480.0) * (n * n * n * n * n);
        double b4 = (4397.0 / 161280.0) * (n * n * n * n) - (11.0 / 504.0) * (n * n * n * n * n);
        double b5 = (4583.0 / 161280.0) * (n * n *        n * n * n * n * n);
        return std::vector<double>{b0, b1, b2, b3, b4, b5};
    };

    auto delta_array = [](double n) -> std::vector<double> {
        double d0 = std::nan(""); // dummy
        double d1 = 2.0 * n - (2.0 / 3.0) * (n * n) - 2.0 * (n * n * n) + (116.0 / 45.0) * (n * n * n * n) + (26.0 / 45.0) * (n * n * n * n * n) - (2854.0 / 675.0) * (n * n * n * n * n * n);
        double d2 = (7.0 / 3.0) * (n * n) - (8.0 / 5.0) * (n * n * n) - (227.0 / 45.0) * (n * n * n * n) + (2704.0 / 315.0) * (n * n * n * n * n) + (2323.0 / 945.0) * (n * n * n * n * n * n);
        double d3 = (56.0 / 15.0) * (n * n * n) - (136.0 / 35.0) * (n * n * n * n) - (1262.0 / 105.0) * (n * n * n * n * n) + (73814.0 / 2835.0) * (n * n * n * n * n * n);
        double d4 = (4279.0 / 630.0) * (n * n * n * n) - (332.0 / 35.0) * (n * n * n * n * n) - (399572.0 / 14175.0) * (n * n * n * n * n * n);
        double d5 = (4174.0 / 315.0) * (n * n * n * n * n) - (144838.0 / 6237.0) * (n * n * n * n * n * n);
        double d6 = (601676.0 / 22275.0) * (n * n * n * n * n * n);
        return std::vector<double>{d0, d1, d2, d3, d4, d5, d6};
    };

    // 定数 (a, F: 世界測地系-測地基準系1980（GRS80）楕円体)
    double m0 = 0.9999;
    double a = 6378137.0;
    double F = 298.257222101;

    // (1) n, A_i, beta_i, delta_iの計算
    double n = 1.0 / (2 * F - 1);
    std::vector<double> A_arr = A_array(n);
    std::vector<double> beta_arr = beta_array(n);
    std::vector<double> delta_arr = delta_array(n);

    // (2), S, Aの計算
    double A_ = ((m0 * a) / (1.0 + n)) * A_arr[0];
    double S_ = ((m0 * a) / (1.0 + n)) * (A_arr[0] * phi0_rad + std::inner_product(A_arr.begin() + 1, A_arr.end(), std::begin(std::vector<double>{2 * phi0_rad, 4 * phi0_rad, 6 * phi0_rad, 8 * phi0_rad, 10 * phi0_rad}), 0.0));

    // (3) xi, etaの計算
    double xi = (x + S_) / A_;
    double eta = y / A_;

    // (4) xi', eta'の計算
    double xi2 = xi;
    double eta2 = eta;
    for (int i = 1; i <= 5; ++i) {
        xi2 -= beta_arr[i] * std::sin(2 * i * xi) * std::cosh(2 * i * eta);
        eta2 -= beta_arr[i] * std::cos(2 * i * xi) * std::sinh(2 * i * eta);
    }

    // (5) chiの計算
    double chi = std::asin(std::sin(xi2) / std::cosh(eta2)); // [rad]
    double latitude = chi + std::inner_product(delta_arr.begin() + 1, delta_arr.end(), std::begin(std::vector<double>{2 * chi, 4 * chi, 6 * chi, 8 * chi, 10 * chi, 12 * chi}), 0.0); // [rad]

    // (6) 緯度(latitude), 経度(longitude)の計算
    double longitude = lambda0_rad + std::atan(std::sinh(eta2) / std::cos(xi2)); // [rad]

    // ラジアンを度になおしてreturn
    return std::make_pair(latitude * 180.0 / M_PI, longitude * 180.0 / M_PI); // [deg]
}



