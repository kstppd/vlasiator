#include "genericTsPool.h"
#include "matrix.h"
#include "tinyAI.h"
#include <algorithm>
#include <array>
#include <vector>

// #define USE_GPU

typedef double Real;
typedef float Realf;
using namespace NumericMatrix;
using namespace TINYAI;

template <typename T> T* check_ptr(T* ptr) {
   if (ptr) {
      return ptr;
   }
   std::cerr << "INVALID POINTER DETECTED" << std::endl;
   abort();
   return nullptr;
}

void scale_vdf(MatrixView<Real>& vspace, Real sparse) {
   constexpr Real minValue = static_cast<Real>(0.001);
   std::for_each(vspace.begin(), vspace.end(),
                 [sparse](Real& value) { value = std::abs(std::log10(std::max(value, minValue * sparse))); });
}

std::array<Real, 2> normalize_vdf(MatrixView<Real>& vdf) {
   Real min_val = *(check_ptr(std::min_element(vdf.begin(), vdf.end())));
   Real max_val = *(check_ptr(std::max_element(vdf.begin(), vdf.end())));
   Real range = max_val - min_val;

   std::for_each(vdf.begin(), vdf.end(), [min_val, range](Real& value) { value = (value - min_val) / range; });
   return {min_val, max_val};
}

void unnormalize_vdf(HostMatrix<Real>& vdf, std::array<Real, 2> norm) {
   Real max_val = norm[1];
   Real min_val = norm[0];
   Real range = max_val - min_val;

   std::for_each(vdf.begin(), vdf.end(), [min_val, range](Real& value) { value = value * range + min_val; });
}

void unscale_vdf(HostMatrix<Real>& vdf) {
   std::for_each(vdf.begin(), vdf.end(), [](Real& value) { value = std::pow(10.0, -1.0 * value); });
}

void sparsify(HostMatrix<Real>& vdf, Real sparse) {
   return;
   std::for_each(vdf.begin(), vdf.end(), [sparse](Real& x) {
      if (x - sparse <= 0.0) {
         x = 0.0;
      }
   });
}

template <BACKEND HW>
NumericMatrix::Matrix<Real, HW> add_fourier_features(const MatrixView<Real>& vcoords, std::size_t order,
                                                     std::vector<Real>& harmonics, GENERIC_TS_POOL::MemPool* p) {
   if (harmonics.empty()) {
      harmonics.resize(order);
      std::random_device rd;
      std::mt19937 gen(rd());
      std::normal_distribution<Real> dist(0, 24);
      std::generate(harmonics.begin(), harmonics.end(), [&]() { return std::abs(dist(gen)); });
   }
   const size_t totalDims = 3 + order * 6;

   NumericMatrix::HostMatrix<Real> host_encoded_vspace(vcoords.nrows(), totalDims);
   for (std::size_t i = 0; i < vcoords.nrows(); ++i) {
      Real vx = vcoords(i, 0) - 0.5;
      Real vy = vcoords(i, 1) - 0.5;
      Real vz = vcoords(i, 2) - 0.5;
      assert(vx >= -0.5 && vx <= 0.5);
      assert(vy >= -0.5 && vx <= 0.5);
      assert(vz >= -0.5 && vx <= 0.5);
      host_encoded_vspace(i, 0) = vx;
      host_encoded_vspace(i, 1) = vy;
      host_encoded_vspace(i, 2) = vz;
      for (std::size_t f = 0; f < order; ++f) {
         host_encoded_vspace(i, f * 6 + 3) = std::sin(harmonics[f] * 2.0 * M_PI * vx);
         host_encoded_vspace(i, f * 6 + 4) = std::sin(harmonics[f] * 2.0 * M_PI * vy);
         host_encoded_vspace(i, f * 6 + 5) = std::sin(harmonics[f] * 2.0 * M_PI * vz);
         host_encoded_vspace(i, f * 6 + 6) = std::cos(harmonics[f] * 2.0 * M_PI * vx);
         host_encoded_vspace(i, f * 6 + 7) = std::cos(harmonics[f] * 2.0 * M_PI * vy);
         host_encoded_vspace(i, f * 6 + 8) = std::cos(harmonics[f] * 2.0 * M_PI * vz);
      }
   }

   NumericMatrix::Matrix<Real, HW> device_encoded_vspace(host_encoded_vspace.nrows(), host_encoded_vspace.ncols(), p);
   NumericMatrix::get_from_host(device_encoded_vspace, host_encoded_vspace);
   return device_encoded_vspace;
}

std::size_t compress_and_reconstruct_vdf(const MatrixView<Real>& vcoords, const MatrixView<Real>& vspace,
                                  MatrixView<Real>& inference_coords, std::size_t fourier_order, std::size_t max_epochs,
                                  std::vector<int>& arch, Real tolerance, HostMatrix<Real>& reconstructed_vdf) {

   constexpr size_t MEMPOOL_BYTES = 2ul * 1024ul * 1024ul * 1024ul;
#ifdef USE_GPU
   constexpr auto HW = BACKEND::DEVICE;
   void* mem;
   cudaMalloc(&mem, MEMPOOL_BYTES);
   std::cout << mem << std::endl;
#else
   constexpr auto HW = BACKEND::HOST;
   void* mem = (void*)malloc(MEMPOOL_BYTES);
#endif
   GENERIC_TS_POOL::MemPool p(mem, MEMPOOL_BYTES);
   std::size_t network_size=0;
   {
      // DO NOT DELETE
      std::vector<Real> harmonics = {2.8370530569956656, 0.06317259286784394, 2.87033597001838,
                                     5.270843933553623,  1.7121147529026062,  0.4102272506250313};
      // std::vector<Real> harmonics;
      NumericMatrix::Matrix<Real, HW> vcoords_train = add_fourier_features<HW>(vcoords, fourier_order, harmonics, &p);
      NumericMatrix::Matrix<Real, HW> vspace_train(vspace.nrows(), vspace.ncols(), &p);
      NumericMatrix::Matrix<Real, HW> vcoords_inference =
          add_fourier_features<HW>(inference_coords, fourier_order, harmonics, &p);
      NumericMatrix::Matrix<Real, HW> vspace_inference(inference_coords.nrows(), 1, &p);
      // Actually read in the vspace for training
      if constexpr (HW == BACKEND::HOST) {
         vspace_train.copy_to_host_from_host_view(vspace);
      } else {
         vspace_train.copy_to_device_from_host_view(vspace);
      }

      constexpr size_t BATCHSIZE = 32;
      NeuralNetwork<Real, HW> nn(arch, &p, vcoords_train, vspace_train, BATCHSIZE);
      network_size=nn.get_network_size();

      for (std::size_t i = 0; i < max_epochs; i++) {
         auto l = nn.train(BATCHSIZE, 1.0e-4);
         // if (i % 1 == 0) {
         //    printf("Loss at epoch %zu: %f\n", i, l);
         // }
         if (l < tolerance) {
            break;
         }
      }
      tinyAI_gpuDeviceSynchronize();
      nn.evaluate(vcoords_inference, vspace_inference);
      vspace_inference.export_to_host(reconstructed_vdf);
   }
   // p.defrag();
   // printf("Pool HWM = %f\n", p.memory_hwm());
   // p.stats();
#ifdef USE_GPU
   cudaFree(mem);
#else
   free(mem);
#endif
   return network_size;
}

extern "C" {
Real compress_and_reconstruct_vdf_2(std::array<Real, 3>* vcoords_ptr, Realf* vspace_ptr, std::size_t size,
                                    std::array<Real, 3>* inference_vcoords_ptr, Realf* new_vspace_ptr,
                                    std::size_t inference_size, std::size_t max_epochs, std::size_t fourier_order,
                                    size_t* hidden_layers_ptr, size_t n_hidden_layers, Real sparsity, Real tol,
                                    Real* weights_ptr, std::size_t weight_size, bool use_input_weights) {

   std::vector<Real> vdf;
   vdf.reserve(size);
   for (std::size_t i = 0; i < size; ++i) {
      vdf.push_back(static_cast<Real>(vspace_ptr[i]));
   }

   const std::size_t vdf_size=vdf.size()*sizeof(Real);

   std::vector<int> arch;
   arch.reserve(n_hidden_layers + 1);
   for (size_t i = 0; i < n_hidden_layers; ++i) {
      arch.push_back(static_cast<int>(hidden_layers_ptr[i]));
   }
   arch.push_back(1);

   MatrixView<Real> vcoords = get_view_from_raw(&(vcoords_ptr[0][0]), size, 3);
   MatrixView<Real> inference_coords = get_view_from_raw(&(inference_vcoords_ptr[0][0]), inference_size, 3);
   MatrixView<Real> vspace = get_view_from_raw(vdf.data(), vdf.size(), 1);
   HostMatrix<Real> vspace_inference_host(inference_coords.nrows(), 1);

   // Scale and normalize
   scale_vdf(vspace, sparsity);
   std::array<Real, 2> norm = normalize_vdf(vspace);

   // Reconstruct
   const std::size_t bytes_used=compress_and_reconstruct_vdf(vcoords, vspace, inference_coords, fourier_order, max_epochs, arch, tol,
                                vspace_inference_host);

   // Undo scalings
   unnormalize_vdf(vspace_inference_host, norm);
   unscale_vdf(vspace_inference_host);
   sparsify(vspace_inference_host, sparsity);

   // Copy back
   for (std::size_t i = 0; i < vspace_inference_host.nrows(); ++i) {
      new_vspace_ptr[i] = static_cast<Realf>(vspace_inference_host(i, 0));
   }
   return static_cast<float>(vdf_size)/static_cast<float>(bytes_used);
}
}
