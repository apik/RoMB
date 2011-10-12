#include "uf.h"
//#include "utils.h"
#include <boost/lexical_cast.hpp>
#include <boost/fusion/container/generation/make_map.hpp>
using std::string;
using std::cout;
using std::endl;


const symbol & get_symbol(const string & s)
{
  static std::map<string, symbol> directory;
  std::map<string, symbol>::iterator i = directory.find(s);
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(make_pair(s, symbol(s))).first->second;
}


UFXmap UF(lst k_lst,lst p_lst,lst subs_lst, unsigned int displacement )
{
  /// constructing expr Xi*Pi
  ex fprop = 0;
  ex sDtmp;
  string str;
  lst xp;
  for(lst::const_iterator it = p_lst.begin();it!=p_lst.end();++it)
    {
      str = "x_" + boost::lexical_cast<string>(displacement+std::distance(p_lst.begin(),it));
      xp.append(get_symbol(str)); // storing list of X identities
      fprop += get_symbol(str) * (*it); 
    }
  fprop = fprop.expand();
  sDtmp = fprop;
  //cout<<fprop<<endl;
  matrix M(k_lst.nops(),k_lst.nops());
  matrix Q(k_lst.nops(),1);//column
  GiNaC::numeric half(1,2);
  for(lst::const_iterator itr = k_lst.begin();itr!=k_lst.end();++itr)
    for(lst::const_iterator itc = itr;itc!=k_lst.end();++itc)
      if(itr == itc)
        {
          M(distance(k_lst.begin(),itr),distance(k_lst.begin(),itc)) = -1*fprop.coeff((*itr)*(*itc),1);
          //cout<<(*itr)*(*itc)<<"M("<<distance(k_lst.begin(),itr)<<","<<distance( k_lst.begin(),itc)<<") coeff "<<fprop.coeff((*itr)*(*itc),1)<<endl;
          sDtmp -= (*itr)*(*itc)*fprop.coeff((*itr)*(*itc),1);
        }
      else
        {
          M(distance(k_lst.begin(),itr),distance( k_lst.begin(),itc)) = -1*half*fprop.coeff((*itr),1).coeff((*itc),1);
          M(distance(k_lst.begin(),itc),distance(k_lst.begin(),itr)) = -1*half*fprop.coeff((*itr),1).coeff((*itc),1);
          //cout<<(*itr)*(*itc)<<"M("<<distance( k_lst.begin(),itr)<<","<<distance( k_lst.begin(),itc)<<") coeff "<<fprop.coeff((*itr),1).coeff((*itc),1)<<endl;
          sDtmp -= (*itr)*(*itc)*fprop.coeff((*itr),1).coeff((*itc),1);
        }
  //cout<<"M: "<<M<<endl;
  sDtmp = sDtmp.expand();
  //cout<<"Expr linear on external momentum: "<<sDtmp.expand()<<endl;

  for(lst::const_iterator itr = k_lst.begin();itr!=k_lst.end();++itr)
    {
      Q(distance( k_lst.begin(),itr),0) = half*sDtmp.coeff((*itr),1);
      sDtmp -= (*itr)*sDtmp.coeff((*itr),1);
    }
  //  cout<<"Q: "<<Q<<endl;
  sDtmp = sDtmp.expand();
  ex minusJ = sDtmp;
   cout<<"-J: "<<minusJ<<endl;
  ex U = M.determinant();
  ex F = expand(M.determinant()*(minusJ+Q.transpose().mul(M.inverse()).mul(Q)(0,0)));
  lst lp;
  F=F.normal();
  cout<<"U= "<<U<<endl<<"F= "<<F<<endl;
  //  cout<<"pol_list "<<lp<<endl;
  U=U.subs(subs_lst,subs_options::algebraic);
  F=F.subs(subs_lst,subs_options::algebraic);
  //  cout<<"UF_WORK"<<endl;
  return fusion::make_map<UFX::U,UFX::F,UFX::xlst>(U,F,xp);
}
