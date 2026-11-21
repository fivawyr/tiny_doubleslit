#include "raylib.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <complex>
#include <array>
#include <algorithm>

using namespace std::complex_literals;
using namespace std;

template <typename T>
string to_string_with_precision(const T a_value, const int n = 2) {
    ostringstream out;
    out.precision(n);
    out << fixed << a_value;
    return out.str();
}

Color colormap(double t) {
    t = clamp(t, 0.0, 1.0);
    unsigned char r = static_cast<unsigned char>(255.0 * clamp(1.5 * t - 0.2, 0.0, 1.0));
    unsigned char g = static_cast<unsigned char>(255.0 * clamp(sin(t * 3.14159265), 0.0, 1.0));
    unsigned char b = static_cast<unsigned char>(255.0 * clamp(1.0 - 0.7 * t, 0.0, 1.0));
    return Color{ r, g, b, 255 };
}

Image density_to_image(const mat& density, const mat* potential = nullptr) {
    const int width = static_cast<int>(density.n_cols);
    const int height = static_cast<int>(density.n_rows);

    double max_val = density.max();
    if (max_val <= 1e-15) max_val = 1.0;

    Image img = GenImageColor(width, height, BLACK);
    Color* pixels = static_cast<Color*>(img.data);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            double norm = density(y, x) / max_val;
            Color c = colormap(std::sqrt(norm));

            if (potential && (*potential)(y, x) > 1.0) {
                c = Color{ 230, 230, 240, 255 };
            }
            pixels[y * width + x] = c;
        }
    }
    return img;
}

void view_slice_raylib(const mat& density, const mat* potential = nullptr, const string& title = "Wavepacket") {
    const int tex_w = static_cast<int>(density.n_cols);
    const int tex_h = static_cast<int>(density.n_rows);

    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    InitWindow(800, 800, title.c_str());
    SetTargetFPS(60);

    Image img = density_to_image(density, potential);
    Texture2D tex = LoadTextureFromImage(img);
    UnloadImage(img);
    SetTextureFilter(tex, TEXTURE_FILTER_BILINEAR);

    bool show_potential = (potential != nullptr);

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_P) && potential != nullptr) {
            show_potential = !show_potential;
            Image updated_img = density_to_image(density, show_potential ? potential : nullptr);
            UpdateTexture(tex, updated_img.data);
            UnloadImage(updated_img);
        }

        BeginDrawing();
        ClearBackground(BLACK);

        Rectangle src = { 0.0f, 0.0f, static_cast<float>(tex_w), static_cast<float>(tex_h) };
        Rectangle dst = {
            0.0f,
            0.0f,
            static_cast<float>(GetScreenWidth()),
            static_cast<float>(GetScreenHeight())
        };
        DrawTexturePro(tex, src, dst, Vector2{ 0.0f, 0.0f }, 0.0f, WHITE);

        DrawText("ESC: Beenden | P: Potenzial an/aus", 12, 12, 18, RAYWHITE);
        EndDrawing();
    }

    UnloadTexture(tex);
    CloseWindow();
}

void simulate(string inputfile) {
    ifstream myfile(inputfile);
    if (!myfile.is_open()) {
        cout << "Unable to open the file " << inputfile << endl;
        return;
    }

    string line;
    std::getline(myfile, line); // Header-Zeile überspringen
    const size_t input_vals = 14;
    std::vector<std::array<double, input_vals>> vars;

    while (std::getline(myfile, line)) {
        if (line.empty()) continue;
        double prob, h, deltat, T, xc, sx, px, yc, sy, py, v0, reim, last_slice, slits;
        std::stringstream mysstream(line);
        if (mysstream >> prob >> h >> deltat >> T >> xc >> sx >> px >> yc >> sy >> py >> v0 >> slits >> reim >> last_slice) {
            vars.push_back({{prob, h, deltat, T, xc, sx, px, yc, sy, py, v0, slits, reim, last_slice}});
        }
    }
    myfile.close();

    int width = 10;
    cout << std::setw(width) << "Problem" << std::setw(width) << "h" << std::setw(width) << "deltat"
         << std::setw(width) << "T" << std::setw(width) << "x_c" << std::setw(width) << "sigma_x"
         << std::setw(width) << "p_x" << std::setw(width) << "y_c" << std::setw(width) << "sigma_y"
         << std::setw(width) << "p_y" << std::setw(width) << "v_0" << std::setw(width) << "slits"
         << std::setw(width) << "ReIm" << std::setw(width) << "Last_slice" << endl;

    for (size_t i = 0; i < vars.size(); ++i) {
        for (size_t j = 0; j < input_vals; ++j) {
            cout << std::setw(width) << vars[i][j];
        }
        cout << endl;
    }

    cout << "Number of simulations to run is " << vars.size() << "." << endl;

    for (size_t i = 0; i < vars.size(); ++i) {
        double prob = vars[i][0];
        double h = vars[i][1];
        double deltat = vars[i][2];
        double T = vars[i][3];
        double xc = vars[i][4];
        double sx = vars[i][5];
        double px = vars[i][6];
        double yc = vars[i][7];
        double sy = vars[i][8];
        double py = vars[i][9];
        double v0 = vars[i][10];
        double slits = vars[i][11];
        double reim = vars[i][12];
        double last_slice = vars[i][13];

        cout << "\nRunning problem " << to_string_with_precision(prob)
             << " with parameters: h=" << h << ", delta t=" << deltat << ", T=" << T
             << ", x_c=" << xc << ", s_x=" << sx << ", p_x=" << px
             << ", y_c=" << yc << ", s_y=" << sy << ", p_y=" << py
             << ", v_0=" << v0 << ", slits=" << (int)slits
             << ", ReIm=" << reim << ", last slice=" << last_slice << "." << endl;

        Crank crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, slits);
        cx_cube results_cube;
        cx_mat last_slice_mat;

        if (last_slice == 1) {
            last_slice_mat = crank.run_simulation(1);
        } else {
            results_cube = crank.run_simulation();
        }

        if (reim == 1) {
            results
