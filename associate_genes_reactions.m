function stringa_output = associate_genes_reactions(stringa)

%unfortunately Recon grRules have no parentheses, which means that we need
%to find a smart way to make sure that AND is solved before OR when we
%substitute MIN and MAX respectively. This means that in the final expression, the MINs have to be
%calculated before the MAXs.
%To do so, we substitute the ORs first (which become MAXs), and then the ANDs inside
%the MAXs. This is to ensure that we have an expression that first solves
%the  ANDs (which are internal) and then solves the ORs (which are
%external), thus respecting the common rule that AND is solved before OR
   
    while ( ~isempty(findstr(' or ', stringa)) )   %loops until all the AND and OR are not found because they have been substituted by MIN and MAX
    
        i = 1;
       while ( (strcmp(stringa(i:i+3),' or ')==0) )
            i = i+1;        %while it does not find any 'and' and any 'or', keeps scrolling the array
        end
        
        stringa = substitute(stringa,i);
        %stringa = regexprep(stringa,'//','/'); 
    end
 
    while ( ~isempty(findstr(' and ', stringa) ) )  %loops until all the AND and OR are not found because they have been substituted by MIN and MAX
    
        i = 1;
        
        while ( (strcmp(stringa(i:i+4),' and ')==0)    )
            i = i+1;        %while it does not find any 'and' and any 'or', keeps scrolling the array
        end
        
        stringa = substitute(stringa,i);
      %stringa = regexprep(stringa,'//','/'); %eliminates duplicate /
    end
    
    if (isempty(findstr('max(', stringa)) && isempty(findstr('min(', stringa)) && ~strcmp(stringa,''))
        stringa_output = stringa;  
    else
        stringa_output = stringa; %if stringa is empty or is a combination of max/min of genes, we leave it as it is
    end
end
    
    
    



    function stringa = substitute(stringa,i)
         
        
        i = i+1;
        %a questo punto i e' posizionato sulla a di 'and' o sulla o di 'or' (prima era posizionato sullo spazio vuoto che precedeva and oppure or), e lo capiamo con questo if
        if (stringa(i)=='a')
            trovato_and = 1;
        else
            trovato_and = 0 ;
        end
            
        
        parentesi_trovate=0;
        
        j = i;
        while (  (strcmp(stringa(j),'(')==0) || (parentesi_trovate~=-1) ) && (strcmp(stringa(j),',')==0 || (parentesi_trovate~=0) ) && (j>1)   %si ferma se trova una parentesi ( che non e' intermedia. E' necessario questo altrimenti si fermerebbe alla prima parentesi ( trovata, anche se questa non e' quella relativa all'AND, ma appartiene ad esempio a un "sotto-OR" che costituisce il ramo sinistro dell'AND, cioe' nel caso ((C or D) and (B or E))
            j = j-1;
            if (stringa(j)==')')
                parentesi_trovate = parentesi_trovate+1; %segnala altre parentesi trovate lungo il percorso
            end
            if (stringa(j)=='(')
                parentesi_trovate = parentesi_trovate-1; %segnala la chiusura delle parentesi trovate lungo il percorso
            end
        end
        if (parentesi_trovate == -1 || strcmp(stringa(j),',')~=0)    %it means that it exited from the WHILE because it actually found a bracket that was not intermediate, or it found a comma. This IF is necessary because it may happen that it exits because it reaches the end of the stringa, so in this case j is already on the first character of the gene, and j=j+1 must not be executed otherwise the first character will not be copied later in the new stringa
            j = j+1;
        end
        %a questo punto j e' posizionato sul carattere iniziale del primo componente dell'and/or
        
        
        
        if (trovato_and == 1)
            k = i+3;
        else
            k = i+2;
        end
       
        parentesi_trovate=0;
       
        while ( (strcmp(stringa(k),')')==0) || ( parentesi_trovate~=-1 )) && (strcmp(stringa(k),',')==0 || (parentesi_trovate~=0) ) && (k<length(stringa))    %si ferma se trova una parentesi ) che non e' intermedia. E' necessario questo altrimenti si fermerebbe alla prima parentesi ) trovata, anche se questa non e' quella relativa all'AND, ma appartiene ad esempio a un "sotto-OR" che costituisce il ramo sinistro dell'AND, cioe' nel caso ((C or D) and (B or E))
            k = k+1;
            if (stringa(k)=='(')
                parentesi_trovate = parentesi_trovate+1; %segnala altre parentesi trovate lungo il percorso
            end
            if (stringa(k)==')')
                parentesi_trovate = parentesi_trovate-1; %segnala la chiusura delle parentesi trovate lungo il percorso
            end
        end
        if parentesi_trovate == -1 ||  strcmp(stringa(k),',')~=0 %this IF is necessary for the same reason as the step "if parentesi_trovate == -1, j=j+1" a few lines before
            k = k-1;
        end
        %a questo punto k e' posizionato sul carattere finale del secondo componente dell'and/or
        
        
        parentesi_trovate=0;
        
        if (trovato_and == 1)
            stringa_new = strrep( stringa, stringa(j:k), strcat('min(',stringa(j:i-1),',',stringa(i+4:k), ')') );
        else
            stringa_new = strrep( stringa, stringa(j:k), strcat('max(',stringa(j:i-1),',',stringa(i+3:k), ')') );
        end
        stringa = stringa_new;
    end


    