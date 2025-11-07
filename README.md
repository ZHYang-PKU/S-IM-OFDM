# S-IM-OFDM: Superposed Index-Modulated OFDM for ISAC

[![IEEE](https://img.shields.io/badge/IEEE-TVT-blue)](https://ieeexplore.ieee.org/document/10559944)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-orange)]()
[![License](https://img.shields.io/badge/License-MIT-green)]()
[![GitHub](https://img.shields.io/badge/Open--Source-Success-success)]()

MATLAB implementation of **Superposed IM-OFDM (S-IM-OFDM)** - an enhanced OFDM waveform for Integrated Sensing and Communications (ISAC).

## 📖 Paper Info

**Title:** Superposed IM-OFDM (S-IM-OFDM): An Enhanced OFDM for Integrated Sensing and Communications  
**Authors:** Zonghui Yang, Shijian Gao, Xiang Cheng, Liuqing Yang  
**Journal:** IEEE Transactions on Vehicular Technology  
**Year:** 2024 | **Volume:** 73 | **Issue:** 10 | **Pages:** 15832-15836
**Citation:** 
```matlab
@ARTICLE{10559944,
  author={Yang, Zonghui and Gao, Shijian and Cheng, Xiang and Yang, Liuqing},
  journal={IEEE Transactions on Vehicular Technology}, 
  title={Superposed IM-OFDM (S-IM-OFDM): An Enhanced OFDM for Integrated Sensing and Communications}, 
  year={2024},
  volume={73},
  number={10},
  pages={15832-15836},
  doi={10.1109/TVT.2024.3412213}
}
```

## 🚀 Key Features

- **Dual-functional waveform** for simultaneous communication and sensing
- **Doppler pre-compensation** using sensed parameters
- **Fusion estimation** for improved sensing accuracy  
- **Flexible power allocation** between sensing and communication
- **Superior performance** in time-varying channels

## 🛠️ Quick Start

### Prerequisites
- MATLAB R2021a or later

### Run Simulations
```matlab
% Communication performance  
BER_waveform.m;

% Sensing performance
RMSE_waveform.m;
