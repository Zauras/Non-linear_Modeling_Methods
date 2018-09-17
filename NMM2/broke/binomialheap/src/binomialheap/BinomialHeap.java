/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package binomialheap;
import java.util.Scanner;
/**
 *
 * @author Tomas
 */
public class BinomialHeap {
public static class BinomialHeapNode
{
    int reiksme, gylis;
    BinomialHeapNode tevas;
    BinomialHeapNode brolis;
    BinomialHeapNode vaikas;
 
    public BinomialHeapNode(int k) {
        reiksme = k;
        gylis = 0;
        tevas = null;
        brolis = null;
        vaikas = null;        
    }
    
    public BinomialHeapNode reverse(BinomialHeapNode sibl) {
            BinomialHeapNode ret;
            if (brolis != null)
                ret = brolis.reverse(this);
            else
                ret = this;
            brolis = sibl;
            return ret;
    }
    
    public BinomialHeapNode findMinNode(){
            BinomialHeapNode x = this, y = this;
            int min = x.reiksme;
 
            while (x != null) {
                if (x.reiksme < min) {
                    y = x;
                    min = x.reiksme;
                }
                x = x.brolis;
            }
 
            return y;
    }
    
    public BinomialHeapNode findANodeWithKey(int value){
            BinomialHeapNode temp = this, node = null;
 
            while (temp != null) 
            {
                if (temp.reiksme == value) 
                {
                    node = temp;
                    break;
                }
                if (temp.vaikas == null)
                    temp = temp.brolis;
                else 
                {
                    node = temp.vaikas.findANodeWithKey(value);
                    if (node == null)
                        temp = temp.brolis;
                    else
                        break;
                }
            }
 
            return node;
    }
    
    public int getSize(){
        return (1 + ((vaikas == null) ? 0 : vaikas.getSize()) + ((brolis == null) ? 0 : brolis.getSize()));
    }
}
 

public static class BinHeap
{
    private BinomialHeapNode Nodes;
    private int size;
 
    public BinHeap()
    {
        Nodes = null;
        size = 0;
    }
    
    public boolean isEmpty()
    {
        return Nodes == null;
    }
    
    public int getSize()
    {
        return size;
    }
 
    public void makeEmpty()
    {
        Nodes = null;
        size = 0;
    }
 
    public void insert(int value) 
    {
        if (value > 0)
        {
            BinomialHeapNode temp = new BinomialHeapNode(value);
            if (Nodes == null) 
            {
                Nodes = temp;
                size = 1;
            } 
            else 
            {
                unionNodes(temp);
                size++;
            }
        }
    }
 
    private void union(BinomialHeapNode binHeap) 
    {
        BinomialHeapNode temp1 = Nodes, temp2 = binHeap;
 
        while ((temp1 != null) && (temp2 != null)) 
        {
            if (temp1.gylis == temp2.gylis) 
            {
                BinomialHeapNode tmp = temp2;
                temp2 = temp2.brolis;
                tmp.brolis = temp1.brolis;
                temp1.brolis = tmp;
                temp1 = tmp.brolis;
            } 
            else 
            {
                if (temp1.gylis < temp2.gylis) 
                {
                    if ((temp1.brolis == null) || (temp1.brolis.gylis > temp2.gylis)) 
                    {
                        BinomialHeapNode tmp = temp2;
                        temp2 = temp2.brolis;
                        tmp.brolis = temp1.brolis;
                        temp1.brolis = tmp;
                        temp1 = tmp.brolis;
                    }
                    else 
                    {
                        temp1 = temp1.brolis;
                    }
                }
                else 
                {
                    BinomialHeapNode tmp = temp1;
                    temp1 = temp2;
                    temp2 = temp2.brolis;
                    temp1.brolis = tmp;
                    if (tmp == Nodes) 
                    {
                        Nodes = temp1;
                    }
                }
            }
        }
        if (temp1 == null) 
        {
            temp1 = Nodes;
            while (temp1.brolis != null) 
            {
                temp1 = temp1.brolis;
            }
            temp1.brolis = temp2;
        } 
    }
 
    private void unionNodes(BinomialHeapNode binHeap) 
    {
        union(binHeap);
 
        BinomialHeapNode prevTemp = null, temp = Nodes, nextTemp = Nodes.brolis;
 
        while (nextTemp != null) 
        {
            if ((temp.gylis != nextTemp.gylis) || ((nextTemp.brolis != null) && (nextTemp.brolis.gylis == temp.gylis))) 
            {                
                prevTemp = temp;
                temp = nextTemp;
            } 
            else
            {
                if (temp.reiksme <= nextTemp.reiksme) 
                {
                    temp.brolis = nextTemp.brolis;
                    nextTemp.tevas = temp;
                    nextTemp.brolis = temp.vaikas;
                    temp.vaikas = nextTemp;
                    temp.gylis++;
                } 
                else 
                {
                    if (prevTemp == null) 
                    {
                        Nodes = nextTemp;
                    }
                    else 
                    {
                        prevTemp.brolis = nextTemp;
                    }
                    temp.tevas = nextTemp;
                    temp.brolis = nextTemp.vaikas;
                    nextTemp.vaikas = temp;
                    nextTemp.gylis++;
                    temp = nextTemp;
                }
            }
            nextTemp = temp.brolis;
        }
    }
    public int findMinimum() 
    {
        return Nodes.findMinNode().reiksme;
    }

    public void delete(int value) 
    {
        if ((Nodes != null) && (Nodes.findANodeWithKey(value) != null)) 
        {
            decreaseKeyValue(value, findMinimum() - 1);
            extractMin();
        }
    }

    public void decreaseKeyValue(int old_value, int new_value) 
    {
        BinomialHeapNode temp = Nodes.findANodeWithKey(old_value);
        if (temp == null)
            return;
        temp.reiksme = new_value;
        BinomialHeapNode tempParent = temp.tevas;
 
        while ((tempParent != null) && (temp.reiksme < tempParent.reiksme)) 
        {
            int z = temp.reiksme;
            temp.reiksme = tempParent.reiksme;
            tempParent.reiksme = z;
 
            temp = tempParent;
            tempParent = tempParent.tevas;
        }
    }

    public int extractMin() 
    {
        if (Nodes == null)
            return -1;
 
        BinomialHeapNode temp = Nodes, prevTemp = null;
        BinomialHeapNode minNode = Nodes.findMinNode();
 
        while (temp.reiksme != minNode.reiksme) 
        {
            prevTemp = temp;
            temp = temp.brolis;
        }
 
        if (prevTemp == null) 
        {
            Nodes = temp.brolis;
        }
        else
        {
            prevTemp.brolis = temp.brolis;
        }
 
        temp = temp.vaikas;
        BinomialHeapNode fakeNode = temp;
 
        while (temp != null) 
        {
            temp.tevas = null;
            temp = temp.brolis;
        }
 
        if ((Nodes == null) && (fakeNode == null))
        {
            size = 0;
        } 
        else
        {
            if ((Nodes == null) && (fakeNode != null)) 
            {
                Nodes = fakeNode.reverse(null);
                size = Nodes.getSize();
            }
            else
            {
                if ((Nodes != null) && (fakeNode == null))
                {
                    size = Nodes.getSize();
                }
                else
                {
                    unionNodes(fakeNode.reverse(null));
                    size = Nodes.getSize();
                }
            }
        }
 
        return minNode.reiksme;
    }
    
}  
 public static void DisplayHeap(BinomialHeapNode pri){
        if (pri!=null) {
        DisplayHeap(pri.vaikas);
        System.out.print(pri.reiksme + " ");
        if (pri.tevas == null) System.out.println();
        DisplayHeap(pri.brolis);
        }
 }
    public static void main(String[] args) {

        /*Scanner scan = new Scanner(System.in);
        BinHeap bh = new BinHeap( );
        char ch;
        do    
        {
            System.out.println("\nVeiksmai su Heap nr. 1:\n");
            System.out.println("1. Iterpti ");
            System.out.println("2. Istrinti ");
            System.out.println("3. Dydis");            
            System.out.println("4. Ar tuscias?");
            System.out.println("5. Isvalyti");
 
            int pasirinko = scan.nextInt();            
            switch (pasirinko)
            {
            case 1 : 
                System.out.println("Iveskite sveika skaiciu, kuri norite iterpti");
                bh.insert( scan.nextInt() ); 
                break;                          
            case 2 :     
                System.out.println("Iveskite elementa, kuri norite pasalinti");            
                bh.delete( scan.nextInt() );   
                break;                         
            case 3 : 
                System.out.println("Heap dydis = "+ bh.getSize());
                break;                                   
            case 4 : 
                System.out.println("Heap tuscias? = "+ bh.isEmpty());
                break; 
            case 5 : 
                bh.makeEmpty();
                System.out.println("Heap isvalytas\n");
                break;         
            default :  
                System.out.println("Toks pasirinkimas negalimas \n ");
                break;   
            }
            DisplayHeap(bh.Nodes);
            System.out.println();
            System.out.println("\n Testi? T/N \n");
            ch = scan.next().charAt(0);                        
        } while (ch == 'T'|| ch == 't'); 
        BinHeap bh1=new BinHeap();
        do    
        {
            System.out.println("\nVeiksmai su Heap nr. 2:\n");
            System.out.println("1. Iterpti ");
            System.out.println("2. Istrinti ");
            System.out.println("3. Dydis");            
            System.out.println("4. Ar tuscias?");
            System.out.println("5. Isvalyti");
 
            int pasirinko = scan.nextInt();           
            switch (pasirinko)
            {
            case 1 : 
                System.out.println("Iveskite sveika skaiciu, kuri norite iterpti");
                bh1.insert( scan.nextInt() ); 
                break;                          
            case 2 :     
                System.out.println("Iveskite elementa, kuri norite pasalinti");            
                bh1.delete( scan.nextInt() );   
                break;                         
            case 3 : 
                System.out.println("Heap dydis = "+ bh1.getSize());
                break;                                   
            case 4 : 
                System.out.println("Heap tuscias? = "+ bh1.isEmpty());
                break; 
            case 5 : 
                bh.makeEmpty();
                System.out.println("Heap isvalytas\n");
                break;         
            default :  
                System.out.println("Toks pasirinkimas negalimas \n ");
                break;   
            }
            DisplayHeap(bh1.Nodes);
            System.out.println();
            System.out.println("\n Testi? T/N \n");
            ch = scan.next().charAt(0);                        
        } while (ch == 'T'|| ch == 't'); 
        
        bh1.unionNodes(bh.Nodes);
        DisplayHeap(bh1.Nodes);
        */
                
        System.out.println("Įterpimas įrašų (Įterpimo efektyvumo palyginimas):");
        
        BinHeap bh1 = new BinHeap( );
        BinHeap bh2 = new BinHeap( );
        BinHeap bh3 = new BinHeap( );
        int n=1500000;
        long lStartTime = System.currentTimeMillis();
        for (int i =1;i<n;i++){
            bh1.insert(i);
        }
        long lEndTime = System.currentTimeMillis();
 	long difference = lEndTime - lStartTime;
        System.out.println("Įterpimo trukmė (vienas BinHeap): " + difference);
        lStartTime = System.currentTimeMillis();
        for (int i =1;i<n/2;i++){
            bh2.insert(i);
        }
        lEndTime = System.currentTimeMillis();
        for (int i =n/2;i<n;i++){
            bh3.insert(i);
        }
        long lEndTime2 = System.currentTimeMillis();
        lEndTime=lEndTime>lEndTime2?lEndTime:lEndTime2;
        difference = lEndTime - lStartTime;
        lStartTime = System.currentTimeMillis();
        bh2.unionNodes(bh3.Nodes);
 	lEndTime = System.currentTimeMillis();
        difference = difference+ lEndTime - lStartTime;
        System.out.println("Įterpimo trukmė (2 BIN HEAP su UNION (per 2 procesorius)): " + difference);
        
    }   
        
        
}
