#include "stdafx.h"
#include "reaction.h"

using namespace std;

Reaction::Reaction(string nome, bool r, bool k, bool m){
	name = nome; 
	rev = r; 
	kin = k;
	mes = m;
	var = 0; 
	subs = 0; 
	prod = 0; 
	fc = 0;
	cs = 0;
	cp = 0;
};

void Reaction::list_enz_complexes(){    //verificar coeficiente estequiometrico
		string pname = name;
		enz_comp.push_back(pname);
		int j = 1;
		for (int i = 0; i < subs; i++){
		//while (i < subs){
			//if (form_complex[i]){
				for (int k = 0; k < abs(metabolites[i].second); k++){
					pname = pname + "." + metabolites[i].first;
					enz_comp.push_back(pname);
					cs++;
				}
				j++;
			//}
		//	i++;
		}
		pname = name;
		//while (i < prod + subs - 1){
		for (int i = subs; i < prod + subs; i++){
			//if (form_complex[i]){
				for (int k = 0; k < abs(metabolites[i].second); k++){
					pname = pname + "." + metabolites[i].first;
					if (k > 0)
						j--;
				}
			//}
		}
		//enz_comp.push_back(pname);
		int len;
		int rep;
		for (int i = prod + subs - 1; i >= subs; i--){
			//if (form_complex[i]){
			//if (i > subs || abs(metabolites[i].second) > 1){
			if (i > subs)
				rep = abs(metabolites[i].second);
			else
				rep = abs(metabolites[i].second) - 1;
			for (int k = 0; k < rep; k++){
				len = pname.length();
				pname = pname.substr(0, len - (metabolites[i].first.length() + 1));
				//if (j <= fc){
				if (j <= prod + subs){
					enz_comp.push_back(pname);
					cp++;
					j++;
				}
			}
			//}
			//}
		}
		for (auto it = enz_comp.begin(); it != enz_comp.end();){    //print para teste
			cout << *it++ << endl;                                         //print para teste
			cout << endl;                                                //print para teste
		}                                                                //print para teste
}

void Reaction::decomp_elem_steps(){
	if (kin){
		list_enz_complexes();
		int i = 0;    //contador enz_comp
		int j = 0;    //contador metabolites
		int l = 0;    //contador metabolito com coeficiente maior que 1
		int n = get_num_elem_steps();
		for (int k = 0; k < n; k = k + 2){
			elem_steps.push_back(name + "_" + to_string(k));
			elem_steps.push_back(name + "_" + to_string(k + 1));
			elem_steps_enz_subs.push_back(pair<int, int>(i, k));
			elem_steps_enz_prod.push_back(pair<int, int>(i, k + 1));
			if (i < cs){
				elem_steps_enz_prod.push_back(pair<int, int>(i + 1, k));
				elem_steps_enz_subs.push_back(pair<int, int>(i + 1, k + 1));
				//while (!form_complex[j])
				//	j++;
				elem_steps_met_subs.push_back(pair<string, int>(metabolites[j].first, k));
				elem_steps_met_prod.push_back(pair<string, int>(metabolites[j].first, k + 1));
				i++; 
				if (abs(metabolites[j].second) > 1 && abs(metabolites[j].second) - 1 > l)
					l++;
				else{
					j++;
					l = 0;
				}
			}
			else{
				if (i == cs)
					j = metabolites.size() - 1;
				elem_steps_enz_prod.push_back(pair<int, int>((k + 1 == n - 1) ? 0 : i + 1, k));
				elem_steps_enz_subs.push_back(pair<int, int>((k + 1 == n - 1) ? 0 : i + 1, k + 1));
				//while (!form_complex[j])
				//	j--;
				elem_steps_met_prod.push_back(pair<string, int>(metabolites[j].first, k));
				elem_steps_met_subs.push_back(pair<string, int>(metabolites[j].first, k + 1));
				i++;
				if (abs(metabolites[j].second) > 1 && abs(metabolites[j].second) - 1 > l)
					l++;
				else{
					j--;
					l = 0;
				}
			}
		}
		// for (auto it = elem_steps.begin(); it != elem_steps.end(); it++){          //print para teste
		// 	cout << *it << endl;                                                   //print para teste
		// }                                                                          //print para teste
		// cout << endl;                                                              //print para teste
		// for (auto it = elem_steps_enz_subs.begin(); it != elem_steps_enz_subs.end(); it++){          //print para teste
		// 	cout << "(" << it->first << ", " << it->second << ") " << enz_comp[it->first] << endl;          //print para teste
		// }                                                                          //print para teste
		// cout << endl;                                                              //print para teste
		// for (auto it = elem_steps_enz_prod.begin(); it != elem_steps_enz_prod.end(); it++){          //print para teste
		// 	cout << "(" << it->first << ", " << it->second << ") " << enz_comp[it->first] << endl;       //print para teste
		// }                                                                          //print para teste
		// cout << endl;                                                              //print para teste
		// for (auto it = elem_steps_met_subs.begin(); it != elem_steps_met_subs.end(); it++){          //print para teste
		// 	cout << "(" << it->first << ", " << it->second << ")" << endl;                           //print para teste
		// }                                                                          //print para teste
		// cout << endl;                                                              //print para teste
		// for (auto it = elem_steps_met_prod.begin(); it != elem_steps_met_prod.end(); it++){          //print para teste
		// 	cout << "(" << it->first << ", " << it->second << ")" << endl;                           //print para teste
		// }                                                                          //print para teste
		// cout << endl;                                                              //print para teste
	}
}


void Reaction::add_metabolite(string met, double coef, bool f){
	metabolites.push_back(pair<string, double>(met, coef));
	form_complex.push_back(f);
	if (f)
		fc++;
	if (coef < 0)
		subs++;
	else
		prod++;
}

//void Reaction::print_reaction(ofstream &file){
//	int size = metabolites.size();
//	int i = 0;
//	for (auto it = metabolites.begin(); it != metabolites.end(); it++){
//		file << "(" << it->second << ") " << it->first.get_name() << " ";
//		if (i < size - 1)
//			file << "+ ";
//		i++;
//	}
//	file << "= 0" << endl;
//}

double Reaction::get_coefficient(string met){
	for (auto it = metabolites.begin(); it != metabolites.end(); it++){
		if (it->first == met){
			return it->second;
		}
	}
	return 0.0;
}