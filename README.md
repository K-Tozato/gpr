# Gaussian Process Regression(ガウス過程回帰)   

 - ガウス過程回帰を用いて補間や応答曲面を構築するためのコード   
 - 5種類のカーネル関数を比較可能．  
 - GPR_functions.R: 必要な関数をまとめたコード
 - GPR.main.R: データの読み込み，パラメータ探索，作図を行うコード

## 入力データ

### 学習用データ (/data/training.txt)  
 - 学習用の入出力データが含まれるファイル  
 - 最後の列に出力データ，それ以外の列に入力データが並ぶデータ


### 検証用データ (/data/test.txt)  
 - 構築したモデルの精度を検証するためのデータ   
 - データ形式は入力データと同じ   
 - (なくてもよい)

### その他
 - methodでパラメータ探索方法を指定する．（MCMCまたは勾配法(grad)）


## 出力データ

### 学習用ケースの再構築の結果 (/res/res_training.txt)    
 - 学習用データの入力値を構築したモデルに代入した時の出力値  
 - 1行目に平均値，2行目に標準偏差

### 検証用データの出力結果 (/res/res_test.txt)    
 - 検証用データの入力値を代入したときの出力値  
 - 1行目に平均値，2行目に標準偏差

### ハイパーパラメータ(/res/hpar.txt)


## 関数の説明
 - f.ker: カーネル関数(knumでカーネル関数を指定)   
 - f.gpr: ガウス過程回帰の関数（評価したい入力値x2に対応する出力値を返す関数）   
 - f.lh : 尤度関数の計算(ver2)     
 - hpar_○○: 勾配法orMCMC(マルコフ連鎖モンテカルロ法)を用いてハイパーパラメータを決定(ver2)      


## プログラムの流れ
 - Functions: 関数の定義   
 - Read data: データの読みこみ   
 - Find hyper-parameter: ハイパーパラメータの設定     
 - Reconstruction: 学習データの入力値を構築したモデルで評価     
 - Result for test data: 検証用データの結果を出力    
 - Compare GPR results with data (Plot): GPRの結果とデータの出力結果の比較   
 - Graph: (入力が1次元or2次元の場合) GPRの結果をグラフ（ばらつきを含む）で表現     



## 出力される図の例    
 - GPRの結果とデータの出力結果の比較   
<img src="./img/compare_output.png" width="50%">

 - 出力値の標準偏差を含んだ形のグラフ     
<img src="./img/output_graph.png" width="50%"> 

<img src="./img/output_graph_3d.png" width="50%">



