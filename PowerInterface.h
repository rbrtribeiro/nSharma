/*
 * PowerInterface.h
 *
 *  Created on: Jul 5, 2017
 *      Author: rr
 */

#ifndef POWERINTERFACE_H_
#define POWERINTERFACE_H_

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include "Pstream.H"

#include <mpi.h>

#include "nSharma_common.h"

typedef int Watt;
#define MPI_WATT MPI_INT

typedef int KHz;

template<typename T>
struct PairLimits {
	T lower;
	T upper;
};

typedef PairLimits<Watt> WattPair;
typedef PairLimits<KHz> KHzPair;

#define POWER_INTERFACE_COMMIT 0

using namespace Foam;

class PowerInterface {
public:

	int nodeID;
	std::string nodeName;
	Watt nodeTDP;
	KHz baseFrequency;
	//int nCores;
	bool nodeMaster;
	float cap;

	PowerInterface(int node, std::string nodeName_, bool nodeMaster_,
			float cap_) {
		nodeID = node;
		nodeName = nodeName_;
		nodeMaster = nodeMaster_;
		cap = cap_;
	}
	virtual ~PowerInterface() {

	}

	virtual Watt getCurrentWatts() =0;
	virtual Watt setWatts(Watt w) =0;
	virtual void getHWWattLimits(WattPair&) =0;
	virtual void getCappedWattLimits(WattPair&) =0;
	virtual std::string getGovernor()=0;
	virtual void reset()=0;
	virtual void printCurrentFullStatus()=0;
	virtual void setPowerCap()=0;
	virtual void setGovernorUserspace()=0;
	virtual KHz getCurrentFrequency() =0;

	KHz WattToFreq(Watt w) {
		return KHz(round((baseFrequency * w) / (nodeTDP * 1.0)));
	}

	Watt FreqToWatt(KHz f) {
		return Watt(round((f * nodeTDP) / (baseFrequency * 1.0)));
	}

	Watt TDP() {
		return nodeTDP;
	}

	bool IsNodeMaster() {
		return nodeMaster;
	}

	std::string toString(int v) {

		std::stringstream ss;
		ss << v;
		return ss.str();
	}

	template<typename Out>
	void split(const std::string &s, char delim, Out result) {
		std::stringstream ss;
		ss.str(s);
		std::string item;
		while (std::getline(ss, item, delim)) {
			*(result++) = item;
		}
	}

	std::vector<std::string> split(const std::string &s, char delim) {
		std::vector < std::string > elems;
		split(s, delim, std::back_inserter(elems));
		return elems;
	}

private:
	virtual void setGovernor(std::string gov)=0;

};

//class SysFsPowerInterface: public PowerInterface {
//public:
//
//	std::string sysFsRoot;
//
//	SysFsPowerInterface(int node, std::string nodeName, std::string sysFsRoot_,
//			bool nodeMaster_, float cap_) :
//			PowerInterface(node, nodeName, nodeMaster_, cap_) {
//
//		sysFsRoot = sysFsRoot_;
//
//	}
//	virtual ~SysFsPowerInterface() {
//
//	}
//
//	int writeToFile(std::string filename, std::string value) {
//
//		std::cout << "writting " + value + " to " + filename << std::endl;
//
//		std::ofstream file;
//		file.open(filename.c_str());
//		if (file.is_open()) {
//			file << value << std::endl;
//			file.close();
//			return 0;
//		} else {
//			std::cout << "FATAL: COULD NOT WRITE TO FILE" << std::endl;
//			return -1;
//		}
//	}
//
//	std::string readFile(std::string filename) {
//		std::ifstream file(filename.c_str());
//		std::string r;
//		getline(file, r);
//		file.close();
//		return r;
//	}
//
//	virtual Watt getCurrentWatts() {
//		//Assuming a bynode power assignment, all the cores should have the same value
//		std::string cpu0 = sysFsRoot + "/cpu0/cpufreq/scaling_cur_freq";
//		KHz f = atoi(readFile(cpu0).c_str());
//		return FreqToWatt(f);
//
//	}
//
//	virtual Watt setWattsCore(Watt w, int core) {
//		std::string cpu = sysFsRoot + "/cpu" + toString(core)
//				+ "/cpufreq/scaling_setspeed";
//		writeToFile(cpu, toString(WattToFreq(w)));
//	}
//
//	virtual Watt setWatts(Watt w) {
//		for (int c = 0; c < nCores; ++c) {
//			setWattsCore(w, c);
//		}
//	}
//
//	virtual void getWattLimits(WattPair& limits) {
//		//Assuming a bynode power assignment, all the cores should have the same value
//		std::string cpu0 = sysFsRoot + "/cpu0/cpufreq/scaling_min_freq";
//		KHz v = atoi(readFile(cpu0).c_str());
//		limits.lower = FreqToWatt(v);
//
//		cpu0 = sysFsRoot + "/cpu0/cpufreq/scaling_max_freq";
//		v = atoi(readFile(cpu0).c_str());
//		limits.upper = FreqToWatt(v);
//	}
//
//	virtual std::string getGovernor() {
//		//Assuming a bynode power assignment, all the cores should have the same value
//		std::string cpu0 = sysFsRoot + "/cpu0/cpufreq/scaling_governor";
//		return readFile(cpu0);
//	}
//
//	virtual void reset() {
//		setGovernor("ondemand");
//	}
//
//	virtual void printCurrentFullStatus() {
//
//	}
//
//	virtual void setPowerCap() {
//
//	}
//
//	virtual void setGovernorUserspace() {
//
//		if (!IsNodeMaster())
//			return;
//
//		setGovernor("userspace");
//	}
//
//private:
//	virtual void setGovernorCore(std::string gov, int core) {
//		std::string cpu = sysFsRoot + "/cpu" + toString(core)
//				+ "/cpufreq/scaling_governor";
//		writeToFile(cpu, gov);
//	}
//
//	virtual void setGovernor(std::string gov) {
//		for (int c = 0; c < nCores; ++c) {
//			setGovernorCore(gov, c);
//
//		}
//	}
//};

class VoidPowerInterface: public PowerInterface {
public:

	VoidPowerInterface(int node, std::string nodeName_, bool nodeMaster_,
			float cap_) :
			PowerInterface(node, nodeName_, nodeMaster_, cap_) {

	}
	virtual ~VoidPowerInterface() {

	}

	virtual Watt getCurrentWatts() {
		return -1;
	}
	virtual Watt setWatts(Watt w) {
		return -1;
	}
	virtual void getHWWattLimits(WattPair& p) {

	}
	virtual void getCappedWattLimits(WattPair& p) {

	}
	virtual std::string getGovernor() {
		return "";
	}

	virtual void reset() {

	}
	virtual void printCurrentFullStatus() {

	}
	virtual void setPowerCap() {

	}
	virtual void setGovernorUserspace() {

	}
	virtual KHz getCurrentFrequency() {
		return -1;
	}

private:
	virtual void setGovernor(std::string gov) {

	}

};

class cpupowerPowerInterface: public PowerInterface {
public:

	std::string sysFsRoot;

	cpupowerPowerInterface(int node, std::string nodeName,
			std::string sysFsRoot_, bool nodeMaster_, float cap_) :
			PowerInterface(node, nodeName, nodeMaster_, cap_) {

		std::cout << "Instancing cpupowerPowerInterface" << std::endl;

		sysFsRoot = sysFsRoot_;

	}
	virtual ~cpupowerPowerInterface() {

		if (!IsNodeMaster())
			return;

		reset();
		printCurrentFullStatus();

	}

	std::string exec(std::string cmd) {

//		std::cout << "nSharma::PowerInterface: Node " << nodeID << "("
//						<< nodeName << ") exec: " << cmd << std::endl;

		char buffer[128];
		std::string result = "";
		FILE* pipe = popen(cmd.c_str(), "r");
		if (pipe) {
			while (!feof(pipe)) {
				if (fgets(buffer, 128, pipe) != NULL)
					result += buffer;
			}
		}

		pclose(pipe);
		return result;
	}

	virtual void printCurrentFullStatus() {

		DLM_PRINT_HEADER_START(Pstream::master())

		if (!IsNodeMaster())
			return;\


		std::cout << exec("cpupower frequency-info -o") << std::endl;
		WattPair p;
		getHWWattLimits(p);
		std::cout << "nSharma::PowerInterface: Node " << nodeID << "("
				<< nodeName << ") Power HW limits: " << p.lower << "W, "
				<< p.upper << "W" << std::endl;
		std::cout << exec("cpupower -c all frequency-info -fm") << std::endl;

		DLM_PRINT_HEADER_END(Pstream::master())

	}

	virtual Watt getCurrentWatts() {

		//Assuming a bynode power assignment, all the cores should have the same value
		std::string cpu0 =
				split(exec("cpupower -c 0 frequency-info -f"), '\n')[1];
		KHz f = atoi(cpu0.c_str());
		return FreqToWatt(f);

	}

	virtual KHz getCurrentFrequency() {

		//Assuming a bynode power assignment, all the cores should have the same value
		std::string cpu0 =
				split(exec("cpupower -c 0 frequency-info -f"), '\n')[1];
		KHz f = atoi(cpu0.c_str());
		return f;

	}

	virtual void getHWWattLimits(WattPair& limits) {

		std::string res = exec("cpupower -c 0 frequency-info -l");

		std::vector < std::string > rlimits = split(split(res, '\n')[1], ' ');

		KHz v = atoi(rlimits[0].c_str());
		limits.lower = FreqToWatt(v);

		v = atoi(rlimits[1].c_str());
		limits.upper = FreqToWatt(v);

	}

	virtual void getCappedWattLimits(WattPair& limits) {

		std::string res = exec("cpupower -c 0 frequency-info -l");

		std::vector < std::string > rlimits = split(split(res, '\n')[1], ' ');

		KHz v = atoi(rlimits[0].c_str());
		limits.lower = FreqToWatt(v);

		v = atoi(rlimits[1].c_str());
		limits.upper = FreqToWatt(v) * cap;

	}

	virtual void getAllowedWattLimits(WattPair& limits) {

		std::vector < std::string > rlimits = split(
				split(exec("cpupower -c 0 frequency-info -p"), '\n')[1], ' ');

		KHz v = atoi(rlimits[0].c_str());
		limits.lower = FreqToWatt(v);

		v = atoi(rlimits[1].c_str());
		limits.upper = FreqToWatt(v);

	}

	virtual std::string getGovernor() {

		return split(split(exec("cpupower -c 0 frequency-info -p"), '\n')[1],
				' ')[2];
	}

	/**
	 * Sets
	 */

	virtual void reset() {

		if (!IsNodeMaster())
			return;

		setGovernor("ondemand");
		setAllowedFreqLimits();
	}

	virtual void setGovernorUserspace() {

		if (!IsNodeMaster())
			return;

		setGovernor("userspace");
		setAllowedFreqLimits();
	}

	virtual void setPowerCap() {

		if (!IsNodeMaster())
			return;

		WattPair limits;
		getCappedWattLimits(limits);

#if(POWER_INTERFACE_COMMIT)

		if (cap < 1.0) {
			setWatts(limits.upper);
		}

		DLM_PRINT_HEADER_START(Pstream::master())

		std::cout << "nSharma::PowerInterface: Node " << nodeID << "("
				<< nodeName << ") system power cap set to " << limits.upper
				<< "W" << std::endl;

		DLM_PRINT_HEADER_END(Pstream::master())
#endif

	}

	virtual Watt setWatts(Watt w) {

		if (!IsNodeMaster())
			return 0;

#if(POWER_INTERFACE_COMMIT)

		setGovernorUserspace();

		exec(
				"sudo cpupower -c all frequency-set -f "
						+ toString(WattToFreq(w)));

		DLM_PRINT_HEADER_START(Pstream::master())

		std::cout << "nSharma::PowerInterface: Node " << nodeID << "("
				<< nodeName << ") Watts set to " << getCurrentWatts()
				<< "W at frequency " << getCurrentFrequency() << "KHz"
				<< std::endl;

		DLM_PRINT_HEADER_END(Pstream::master())
#endif

	}

private:
	virtual void setAllowedFreqLimits() {

		std::string res = exec("cpupower -c 0 frequency-info -l");

		std::vector < std::string > rlimits = split(split(res, '\n')[1], ' ');

		KHz lower = atoi(rlimits[0].c_str());
		KHz upper = atoi(rlimits[1].c_str());

#if(POWER_INTERFACE_COMMIT)

		exec("sudo cpupower -c all frequency-set --min " + toString(lower));
		exec("sudo cpupower -c all frequency-set --max " + toString(upper));

		WattPair p;
		getAllowedWattLimits(p);

		std::cout << "nSharma::PowerInterface: Node " << nodeID << "("
				<< nodeName << ") allowed limits set to " << p.lower << ","
				<< p.upper << "W" << std::endl;

#endif

	}

	virtual void setGovernor(std::string gov) {

		if (!IsNodeMaster())
			return;

#if(POWER_INTERFACE_COMMIT)

		exec("sudo cpupower -c all frequency-set -g " + gov);

		std::cout << "nSharma::PowerInterface: Node " << nodeID << "("
				<< nodeName << ") governor set to " << getGovernor()
				<< std::endl;
#endif

	}
};

class SeARCHPowerInterface: public cpupowerPowerInterface {
public:

	SeARCHPowerInterface(int node, std::string nodeName, std::string sysFsRoot_,
			bool nodeMaster_, float cap_) :
			cpupowerPowerInterface(node, nodeName, sysFsRoot_, nodeMaster_,
					cap_) {

		sysFsRoot = sysFsRoot_;

		if (nodeName.find("641") != std::string::npos) {
			nodeTDP = 95;
			baseFrequency = 2601000;
			//nCores = 32;
		} else if (nodeName.find("421") != std::string::npos) {
			nodeTDP = 80;
			baseFrequency = 2262000;
			//nCores = 16;
		} else if (nodeName.find("662") != std::string::npos) {
			nodeTDP = 115;
			baseFrequency = 2400000;
			//nCores = 48;
		} else if (nodeName.find("erform") != std::string::npos) {
			nodeTDP = 80;
			baseFrequency = 2530000;
			//nCores = 8;
		} else if (nodeName.find("002-1") != std::string::npos) {
			nodeTDP = 215;
			baseFrequency = 1301000;
					//nCores = 8;
		} else {
			std::cout << "FATAL: NODE STATS UNKNOWN" << std::endl;
			nodeTDP = 1;
			baseFrequency = 1;

		}

		reset();

		printCurrentFullStatus();
		setPowerCap();

	}

};

#endif /* POWERINTERFACE_H_ */
