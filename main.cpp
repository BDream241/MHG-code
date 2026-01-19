#include"maintenance.h"
#include <random>

int main() {
	int arg = 5;
	string str = "../data/AM";
	//core分解
	if (arg == 0) {
		BiGraph graph(str);
		cout << "max k-core: " << graph.coreDecompose() << endl;
	}
	//alpha、beta分解
	else if (arg == 1) {
		BiGraph graph(str);
		graph.decompose();
	}
	//query
	else if (arg == 2) {
		BiGraph graph(str);
		auto start = chrono::system_clock::now();
		graph.decompose();
		auto end = chrono::system_clock::now();
		chrono::duration<double> time = end - start;
		cout << "decompose time: " << time.count() << endl;
	
		ASG  myAsg(graph);

		start = chrono::system_clock::now();
		MHG myMhG(myAsg.asg, graph.maxK, graph.maxAlpha, graph.maxBeta, graph.n1, graph.n2);
		end = chrono::system_clock::now();
		time = end - start;
		cout << "MHG time: " << time.count() << endl;


		//for (auto i : myMhG.tf)
			//cout << i.end - i.begin << endl;

		string strq = str + "/query.txt";
		FILE* fq = fopen(strq.c_str(), "r");
		bool isLeft;
		int q, alpha, beta;
		vector<int>aList, bList, nList, isList;
		while (fscanf(fq, "%d%d%d%d", &q, &isLeft, &alpha, &beta) != EOF) {
			aList.push_back(alpha);
			bList.push_back(beta);
			nList.push_back(q);
			isList.push_back(isLeft);
		}
		fclose(fq);
		double sum = 0;

		vector<bool>leftResult(graph.n1 + 1, false);
		vector<bool>rightResult(graph.n2 + 1, false);
		for (int i = 0; i < 5; i++) {
			for (int j = 0; j < 100; j++) {
				start = chrono::system_clock::now();
				leftResult.resize(graph.n1 + 1, false);
				rightResult.resize(graph.n2 + 1, false);
	
				myMhG.QueryViaMHG(aList[i * 100 + j], bList[i * 100 + j], nList[i * 100 + j], isList[i * 100 + j], leftResult, rightResult);
				//verifyCom(graph, leftResult, rightResult, aList[i * 100 + j], bList[i * 100 + j], nList[i * 100 + j], isList[i * 100 + j]);
				leftResult.clear();
				
				end = chrono::system_clock::now();
				rightResult.clear();

				time = end - start;
				//cout << "th time: " << time.count() * 1000 << endl;
				//sum += time.count() * 1000;

			}
		}
		cout << sum << "  " << sum / 500.0 << endl;
	}
	else if (arg == 3) {
		BiGraph graph(str);
		graph.decompose();
		ASG  myAsg(graph);
		string adds = str + "/insert.txt";

		auto start = chrono::system_clock::now();
		//Maintenance mt(adds, myAsg, graph);
		Maintenance mt4(adds, myAsg, graph, 16);
		auto end = chrono::system_clock::now();
		chrono::duration<double> time = end - start;
		cout << "16 time: " << time.count() << endl;

		start = chrono::system_clock::now();
		//Maintenance mt(adds, myAsg, graph);
		Maintenance mt5(adds, myAsg, graph, 32);
		end = chrono::system_clock::now();
		time = end - start;
		cout << "32 time: " << time.count() << endl;

		start = chrono::system_clock::now();
		//Maintenance mt(adds, myAsg, graph);
		Maintenance mt6(adds, myAsg, graph, 64);
		end = chrono::system_clock::now();
		time = end - start;
		cout << "64 time: " << time.count() << endl;

		start = chrono::system_clock::now();
		//Maintenance mt(adds, myAsg, graph);
		Maintenance mt7(adds, myAsg, graph, 128);
		end = chrono::system_clock::now();
		time = end - start;
		cout << "128 time: " << time.count() << endl;

		start = chrono::system_clock::now();
		//Maintenance mt(adds, myAsg, graph);
		Maintenance mt8(adds, myAsg, graph, 256);
		end = chrono::system_clock::now();
		time = end - start;
		cout << "256 time: " << time.count() << endl;

		start = chrono::system_clock::now();
		//Maintenance mt(adds, myAsg, graph);
		Maintenance mt9(adds, myAsg, graph, 512);
		end = chrono::system_clock::now();
		time = end - start;
		cout << "512 time: " << time.count() << endl;

		start = chrono::system_clock::now();
		//Maintenance mt(adds, myAsg, graph);
		Maintenance mt10(adds, myAsg, graph, 1024);
		end = chrono::system_clock::now();
		time = end - start;
		cout << "1024 time: " << time.count() << endl;
	}
	else if (arg == 4) {
		BiGraph graph(str);
		graph.decompose();
		ASG  myAsg(graph);

		/*vector<vector<int>>bk(6);;
		
		for (int i = 0; i < myAsg.uIndex.size(); i++) {
			for (int j = 0; j < myAsg.uIndex[i].size(); j++) {
				int alpha = myAsg.asg[myAsg.uIndex[i][j]].alpha;
				int beta= myAsg.asg[myAsg.uIndex[i][j]].beta;

				if (alpha >= 10 && beta >= 2 && bk[1].size() < 10) {
					bk[1].push_back(i);
					break;
				}
				if (alpha >= 10 && beta >= 4 && bk[2].size() < 10) {

					bk[2].push_back(i);
					break;
				}
				if (alpha >= 10 && beta >= 6 && bk[3].size() < 10) {
					bk[3].push_back(i);
					break;
				}
				if (alpha >= 10 && beta >= 8 && bk[4].size() < 10) {
					bk[4].push_back(i);
					break;
				}
				if (alpha >= 10 && beta >= 10 && bk[5].size() < 10) {
					bk[5].push_back(i);
					break;
				}
				
			}
		}

		string svar = str + "/varBeta.txt";
		FILE* fqq = fopen(svar.c_str(), "w");
		for (int i = 1; i < 6; i++) {
			for (int j = 0; j < bk[i].size(); j++) {
				fprintf(fqq, "%d 1 %d %d\n ",bk[i][j],10, i*2 );
			}
		}
		fclose(fqq);

		cout << bk[5].size() << endl;*/

		MHG myMhG(myAsg.asg, graph.maxK, graph.maxAlpha, graph.maxBeta, graph.n1, graph.n2);




		string strq = str + "/query.txt";
		FILE* fq = fopen(strq.c_str(), "r");
		bool isLeft;
		int q, alpha, beta;
		vector<int>aList, bList, nList, isList;
		while (fscanf(fq, "%d%d%d%d", &q, &isLeft, &alpha, &beta) != EOF) {
			aList.push_back(alpha);
			bList.push_back(beta);
			nList.push_back(q);
			isList.push_back(isLeft);
		}
		fclose(fq);

		vector<bool>leftResult(graph.n1 + 1, false);
		vector<bool>rightResult(graph.n2 + 1, false);

		auto start = chrono::system_clock::now();
		auto end = chrono::system_clock::now();
		chrono::duration<double> time = end - start;

		for (int i = 0; i < 1; i++) {
			cout << i << endl;
			
			for (int j = 0; j < 10; j++) {
				leftResult.resize(graph.n1 + 1, false);
				rightResult.resize(graph.n2 + 1, false);
				start = chrono::system_clock::now();
				myMhG.QueryViaMHG(aList[i * 10 + j], bList[i * 10 + j], nList[i * 10 + j], isList[i * 10 + j], leftResult, rightResult);
				end = chrono::system_clock::now();
				time += end - start;
				//verifyCom(graph, leftResult, rightResult, aList[i * 100 + j], bList[i * 100 + j], nList[i * 100 + j], isList[i * 100 + j]);
				leftResult.clear();
				rightResult.clear();
			}
			
			cout << "query time: " << time.count() * 100 << endl;
		}
	}

	else if (arg == 5) {
		BiGraph graph(str);
		graph.decompose();
		ASG  myAsg(graph);
		//MHG myMhG(myAsg.asg, graph.maxK, graph.maxAlpha, graph.maxBeta, graph.n1, graph.n2);
		for (auto i : myAsg.asg)
			if (i.alpha == 5 && i.beta == 48 )
				cout << i.left.size() << endl;

		vector<bool>leftResult(graph.n1 + 1, false);
		vector<bool>rightResult(graph.n2 + 1, false);
		leftResult.resize(graph.n1 + 1, false);
		rightResult.resize(graph.n2 + 1, false);
		myAsg.Query(5, 45, 22664, 1, leftResult, rightResult);
		verifyCom(graph, leftResult, rightResult, 5, 45, 22664, 1);
		leftResult.clear();
		rightResult.clear();
	}

	return 0;
}