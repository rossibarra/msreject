#include <Sequence/SimData.hpp>
#include <cstdio>

namespace Sequence {
	class SimDataTBS : public Sequence::SimData    {
   	private:
      	unsigned totsam;
      	std::string __tbsline;

    	public:
   		std::string tbsline () const { return __tbsline; }
      	virtual int fromfile( FILE * openedfile );
      	SimDataTBS();
 	};
}
