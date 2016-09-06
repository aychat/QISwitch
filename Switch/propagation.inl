#ifndef PROPAGATION_INL
#define PROPAGATION_INL

template<class Function>
void runge_kutta_step(QI_Observable &x, double dt, Function advance)
{
	class Tmp
	{
	public:
		std::vector<QI_Observable*> tmp;
		Tmp(QI_Model& m)
		{
			tmp.resize(4);
			std::for_each(tmp.begin(), tmp.end(), [&](QI_Observable*& o) {o = new QI_Observable(&m); });
		}
		~Tmp() { std::for_each(tmp.begin(), tmp.end(), [&](QI_Observable*& o) {delete o; }); };
		QI_Observable* operator[](size_t i) { return tmp[i]; };
	};
	static Tmp tmp(*x.model);
	//calculating k1
	*tmp[3] = x;
	tmp[2]->SetZero();
	advance((tmp[2]), &x, 1, tmp[0]);
	x.data += dt / 6 * tmp[2]->data;
	tmp[1]->data = tmp[3]->data + dt / 2 * tmp[2]->data;
	//calculating k2
	tmp[2]->SetZero();
	advance((tmp[2]), tmp[1], 1, tmp[0]);
	x.data += dt / 3 * tmp[2]->data;
	tmp[1]->data = tmp[3]->data + dt / 2 * tmp[2]->data;
	//calculating k3
	tmp[2]->SetZero();
	advance((tmp[2]), tmp[1], 1, tmp[0]);
	x.data += dt / 3 * tmp[2]->data;
	tmp[1]->data = tmp[3]->data + dt * tmp[2]->data;
	//calculating k4
	advance(&x, tmp[1], dt / 6, tmp[0]);
}

template<class Function>
double MakeStep(QI_Observable &x, double delta_t, Function advance)
{
	static QI_Observable tmp2(x.model);
	tmp2.data = x.data;
	runge_kutta_step(x, delta_t / 2,advance);
	runge_kutta_step(x, delta_t / 2,advance);
	runge_kutta_step(tmp2, delta_t,advance);
	tmp2 -= x;
	tmp2 *= (1. / 15);
	x += tmp2;
	return tmp2.data.cwiseAbs().maxCoeff();
};



#endif