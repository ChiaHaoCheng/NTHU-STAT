# 專案目的  

與聯華電子公司合作，根據聯華電子-智慧製造部提供的虛擬量測(Virtual Metrology, VM)資料集，開發晶圓規格(SPEC)之預測或描述模型，藉由 RMSE 指標比較候選模型的能力。  
其中資料集包含兩種Data1和Data2，代表不同的晶圓產品。  

![image](https://github.com/ChiaHaoCheng/NTHU-STAT/blob/eaaf5f1334123e6b8560ded97e3dbb3bf1f8e1b7/Practicing%20Statistics/%E8%81%AF%E8%8F%AF%E9%9B%BB%E5%AD%90/%E8%99%9B%E6%93%AC%E9%87%8F%E6%B8%AC%E6%A6%82%E5%BF%B5%E5%9C%96.png)

## Data1 

針對Data1三種晶圓產品的歷史VM資料，建構預測模型，其中我們嘗試以下模型：  
1. 多項線性迴歸 Polynomial regression  
2. 廣義相加模型 Generalized Additive Model  
3. 集成學習 ensemble learning ( Boosting、Bagging )
![image](https://github.com/ChiaHaoCheng/NTHU-STAT/blob/862ff154518be7ef714727aba42c6cb1b0f2f30a/Practicing%20Statistics/%E8%81%AF%E8%8F%AF%E9%9B%BB%E5%AD%90/data1_%E5%B0%81%E9%9D%A2.png)
![image](https://github.com/ChiaHaoCheng/NTHU-STAT/blob/862ff154518be7ef714727aba42c6cb1b0f2f30a/Practicing%20Statistics/%E8%81%AF%E8%8F%AF%E9%9B%BB%E5%AD%90/data1_B_result.png)

## Data2 

針對Data2一種晶圓產品的即時VM資料(具備Large-Scale特性)，建構預測模型，分析流程如下：  
1. 對input VM variables進行資料清洗  
2. 執行多重檢定選出可用的VM變數
3. 建立多項式迴歸、集成學習模型  
4. 嘗試使用Lazy Learning: Gaussian Process以提升ensemble learning模型的預測能力  
5. 比較數個模型優劣  
![image](https://github.com/ChiaHaoCheng/NTHU-STAT/blob/3b0024e97099435f0252897e18d9b0bc4c684700/Practicing%20Statistics/%E8%81%AF%E8%8F%AF%E9%9B%BB%E5%AD%90/data2_%E5%B0%81%E9%9D%A2.png)
![image](https://github.com/ChiaHaoCheng/NTHU-STAT/blob/3b0024e97099435f0252897e18d9b0bc4c684700/Practicing%20Statistics/%E8%81%AF%E8%8F%AF%E9%9B%BB%E5%AD%90/data2_flowchart.png)
![image](https://github.com/ChiaHaoCheng/NTHU-STAT/blob/3b0024e97099435f0252897e18d9b0bc4c684700/Practicing%20Statistics/%E8%81%AF%E8%8F%AF%E9%9B%BB%E5%AD%90/data2_result.png)
