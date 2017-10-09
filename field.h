#pragma once

class field
{
public:
	virtual ~field(void);
	 field(void);
	virtual void calc_field()=0;
};
